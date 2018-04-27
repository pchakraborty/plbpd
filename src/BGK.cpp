#include <memory>
#include <iostream>
#include "BGK.hpp"
#include "tbb/tbb.h"
#include <array>
#include <immintrin.h>

BGK::BGK(const LBModel *lbmodel, const Domain *domain):
    LBDynamics(), _lbmodel(lbmodel), _domain(domain) {
    _tau = 3.0f*_domain->getFluidViscosity() + 0.5f;
    _omega = 1./_tau;
}

BGK::~BGK(){}

void BGK::initialize(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    std::vector<float> nlocal(kdim);
    // at all (domain+boundary+buffer) nodes
    for (auto zl=0; zl<zdim+2; ++zl){
        for (auto yl=0; yl<ydim+2; ++yl){
            for (auto xl=0; xl<xdim+2; ++xl){
                auto rholocal = lattice.rho->at(zl,yl,xl);
                std::array<float, 3> ulocal;
                for (auto i=0; i<3; ++i)
                    ulocal[i] = lattice.u->at(zl,yl,xl,i);
                _getEqlbDist(rholocal, ulocal, nlocal);
                for (auto k=0; k<kdim; ++k)
                    lattice.n->at(zl,yl,xl,k) = nlocal[k];
            }
        }
    }
}

void BGK::_getEqlbDist(const float rholocal, const std::array<float, 3> &ulocal, std::vector<float> &nlocal){
    const auto kdim = _lbmodel->getNumberOfDirections();
    auto usq = ulocal[0]*ulocal[0] + ulocal[1]*ulocal[1] + ulocal[2]*ulocal[2];
    auto w = _lbmodel->getDirectionalWeights();
    auto c = _lbmodel->getLatticeVelocities();
    for (auto k=0; k<kdim; ++k){
        auto k3 = k*3;
        auto cu = c[k3+0]*ulocal[0] + c[k3+1]*ulocal[1] + c[k3+2]*ulocal[2];
        nlocal[k] = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq); // neq
    }
}

void BGK::collideAndStream(Lattice &lattice){
    _serialCollideAndStream(lattice);
}

// This is the workhorse
void BGK::_collideAndStreamOnPlane_ref(size_t zl, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeights();

    const array3f *rho = lattice.rho;
    const array4f *u = lattice.u;
    array4f *n = lattice.n;
    array4f *ntmp = lattice.ntmp;
    
    std::array<float, 3> extForce = {0.0f, 0.0f, 0.0f}; // TODO: Get extForce from Domain
    
    for (auto yl=1; yl<ydim+1; ++yl){
        for (auto xl=1; xl<xdim+1; ++xl){
            auto rholocal = rho->at(zl,yl,xl);
            auto rhoinv = 1.0f/rholocal;
            std::array<float, 3> ueq;
            for (auto i=0; i<3; ++i)
                ueq[i] = u->at(zl,yl,xl,i) + extForce[i]*_tau*rhoinv;
            auto usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
            for (auto k=0; k<kdim; ++k){
                auto k3 = k*3;
                auto cu = c[k3+0]*ueq[0] + c[k3+1]*ueq[1] + c[k3+2]*ueq[2];
                auto neq = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                ntmp->at(zl+c[k3+2],yl+c[k3+1],xl+c[k3+0],k) =
                    (1.0f-_omega)*n->at(zl,yl,xl,k) + _omega*neq;
            }
        }
    }
}

// This is the workhorse
void BGK::_collideAndStreamOnPlane(size_t zl, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    const auto kdimPadded = 32;
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w_orig = _lbmodel->getDirectionalWeights();
    std::vector<float> w(w_orig);
    for (auto i=0; i<kdimPadded-kdim; ++i)
        w.push_back(0.0f);

    const array3f *rho = lattice.rho;
    const array4f *u = lattice.u;
    array4f *n = lattice.n;
    array4f *ntmp = lattice.ntmp;
    
    auto extForce = std::array<float, 3>(); // TODO: Get extForce from Domain
    auto ueq = std::array<float, 3>(); // value-initialized to zero

    std::vector<float> cu(kdimPadded, 0.0f), nprime(kdimPadded, 0.0f);
    // const auto alignment = 64; // 64 Byte alignment
    // auto cu = static_cast<float*>(_mm_malloc(kdimPadded*sizeof(float), alignment));
    // auto nprime = static_cast<float*>(_mm_malloc(kdimPadded*sizeof(float), alignment));
    // for (auto i=0; i<kdimPadded; ++i){
    //     cu[i] = 0.0f;
    //     nprime[i] = 0.0f;
    // }
    
    // AVX2 variables
    __m256 _one = _mm256_set1_ps(1.0f);
    __m256 _three = _mm256_set1_ps(3.0f);
    __m256 _fourPointFive = _mm256_set1_ps(4.5f);
    __m256 _minusOnePointFive = _mm256_set1_ps(-1.5f);
    __m256 __omega = _mm256_set1_ps(_omega);
    __m256 _oneMinusOmega = _mm256_set1_ps(1.0f-_omega);
    
    for (auto yl=1; yl<ydim+1; ++yl){
        for (auto xl=1; xl<xdim+1; ++xl){
            auto rholocal = rho->at(zl,yl,xl);
            auto tau_rhoinv = _tau/rholocal;
            for (auto i=0; i<3; ++i)
                ueq[i] = u->at(zl,yl,xl,i) + extForce[i]*tau_rhoinv;
            auto usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
            for (auto k=0; k<kdim; ++k){
                auto k3 = k*3;
                cu[k] = c[k3+0]*ueq[0] + c[k3+1]*ueq[1] + c[k3+2]*ueq[2];
            }
            auto _rholocal = _mm256_set1_ps(rholocal);
            auto _kusq = _mm256_set1_ps(-1.5f*usq); // -1.5*usq
            // Collision - compute nprime
            for (auto k=0; k<kdimPadded/8; ++k){
                auto k8 = k*8;
                auto _w = _mm256_loadu_ps(&w[k8]);
                auto _cu = _mm256_loadu_ps(&cu[k8]);
                auto _nk = _mm256_loadu_ps(n->get(zl,yl,xl,k8)); // access overflow, but harmless
                auto _cusq = _mm256_mul_ps(_cu, _cu);
                // neq = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                auto _neq = _mm256_fmadd_ps(_three, _cu, _one); // 3*cu(k) + 1
                _neq = _mm256_fmadd_ps(_fourPointFive, _cusq, _neq); // 4.5*cusq(k) + 3*cu(k) + 1
                _neq = _mm256_add_ps(_neq, _kusq); // -1.5*usq + 4.5*cusq(k) + 3*cu(k) + 1
                _neq = _mm256_mul_ps(_neq, _rholocal);
                _neq = _mm256_mul_ps(_neq, _w);
                // nprime(k) = (1.0-_omega)*n(zl,yl,xl,k) + _omega*neq(k);
                auto _nprime = _mm256_fmadd_ps(_oneMinusOmega, _nk, _mm256_mul_ps(__omega, _neq));
                _mm256_storeu_ps(&nprime[k8], _nprime);
            }
            // Streaming
            for (auto k=0; k<kdim; ++k){ // NOTE: 0<=k<=kdim (NOT kdimPadded)
                auto k3 = k*3;
                ntmp->at(zl+c[k3+2],yl+c[k3+1],xl+c[k3+0],k) = nprime[k];
            }
        }
    }

    // _mm_free(cu);
    // _mm_free(nprime);
}

void BGK::_serialCollideAndStream(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    for (auto zl=1; zl<zdim+1; ++zl){
        _collideAndStreamOnPlane(zl, lattice);
    }
    // swap n and ntmp
    array4f *tmp = lattice.n;
    lattice.n = lattice.ntmp;
    lattice.ntmp = tmp;
}

void BGK::_parallelCollideAndStream(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    tbb::parallel_for(size_t(1), zdim+1, [this, &lattice] (size_t zl){
            _collideAndStreamOnPlane(zl, lattice);
     });
    // swap n and ntmp
    array4f *tmp = lattice.n;
    lattice.n = lattice.ntmp;
    lattice.ntmp = tmp;
}

float BGK::getAvgFluidDensity(){
    return 0.0;
}

void BGK::calcMoments(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeights();

    const array4f *n = lattice.n;
    array3f *rho = lattice.rho;
    array4f *u = lattice.u;

    for (auto zl=1; zl<zdim+1; ++zl){
        for (auto yl=1; yl<ydim+1; ++yl){
            for (auto xl=1; xl<xdim+1; ++xl){
                float rholocal = 0.0f;
                std::array<float, 3> ulocal = {0.0f, 0.0f, 0.0f};
                for (auto k=0; k<kdim; ++k){
                    auto nk = n->at(zl,yl,xl,k);
                    rholocal += nk; // \sum_k n[k]
                    auto k3 = k*3;
                    for (auto i=0; i<3; ++i)
                        ulocal[i] += nk*c[k3+i]; // \sum_k n[k]*c[k][i], i=0,1,2
                }
                rho->at(zl,yl,xl) = rholocal;
                // std::cout<<"rho("<<zl<<", "<<yl<<", "<<xl<<"): "<<rholocal<<std::endl;
                auto rhoinv = 1.0f/rholocal;
                for (auto i=0; i<3; ++i)
                    u->at(zl,yl,xl,i) = ulocal[i]*rhoinv;
            }
        }
    }
}

void BGK::_printInfoForDebugging(){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeights();
    
    std::cout<<"xdim: "<<xdim<<", ydim:  "<<ydim<<", zdim: "<<zdim<<std::endl;
    std::cout<<"kdim: "<<kdim<<std::endl;
    std::cout<<"Lattice velocities: ";
    for(auto ic=begin(c); ic!=end(c); ++ic)
        std::cout<<*ic<<" ";
    std::cout<<std::endl;
    std::cout<<"Directional weights: ";
    for(auto iw=begin(w); iw!=end(w); ++iw)
        std::cout<<*iw<<" ";
    std::cout<<std::endl;
}
