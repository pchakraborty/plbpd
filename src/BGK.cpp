#include <memory>
#include <iostream>
#include "BGK.hpp"
#include "tbb/tbb.h"
#include <array>
#include <immintrin.h>
#include <mm_malloc.h>

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
        auto ck = &c[k*3];
        auto cu = ck[0]*ulocal[0] + ck[1]*ulocal[1] + ck[2]*ulocal[2];
        nlocal[k] = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq); // neq
    }
}

void BGK::collideAndStream(Lattice &lattice){
    _serialCollideAndStream(lattice);
}

// This is the workhorse
void BGK::_collideAndStreamOnPlane(size_t zl, Lattice &lattice){
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
                auto ck = &c[k*3];
                auto cu = ck[0]*ueq[0] + ck[1]*ueq[1] + ck[2]*ueq[2];
                auto neq = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                ntmp->at(zl+ck[2],yl+ck[1],xl+ck[0],k) =
                    (1.0f-_omega)*n->at(zl,yl,xl,k) + _omega*neq;
            }
        }
    }
}

// This is the workhorse
void BGK::_collideAndStreamOnPlaneAvx2(size_t zl, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    const auto kdimAvx = _lbmodel->getNumberOfDirectionsAvx();
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeightsAvx();

    // Shorthands - also helps compiler optimize
    const array3f * __restrict__ rho = lattice.rho;
    const array4f * __restrict__ u = lattice.u;
    array4f * __restrict__ n = lattice.n;
    array4f * __restrict__ ntmp = lattice.ntmp;
    
    // External force. TODO: Get extForce from Domain
    auto extForce = std::array<float, 3>();

    // Local variables
    auto ueq = std::array<float, 3>(); // value-initialized to zero
    auto cu = static_cast<float*>(_mm_malloc(kdimAvx*sizeof(float), 64));
    auto nprime = static_cast<float*>(_mm_malloc(kdimAvx*sizeof(float), 64));
    for (auto k=0; k<kdimAvx; ++k){
        cu[k] = 0.0f;
        nprime[k] = 0.0f;
    }
    
    // SIMD variables
    __m256 _one = _mm256_set1_ps(1.0f);
    __m256 _three = _mm256_set1_ps(3.0f);
    __m256 _fourPointFive = _mm256_set1_ps(4.5f);
    __m256 _minusOnePointFive = _mm256_set1_ps(-1.5f);
    __m256 __omega = _mm256_set1_ps(_omega);
    __m256 _oneMinusOmega = _mm256_set1_ps(1.0f-_omega);
    
    for (auto yl=1; yl<ydim+1; ++yl){
        for (auto xl=1; xl<xdim+1; ++xl){

            // ueq and usq
            auto rholocal = rho->at(zl,yl,xl);
            auto tau_rhoinv = _tau/rholocal;
            auto ulocal = u->get(zl,yl,xl,0);
            for (auto i=0; i<3; ++i)
                ueq[i] = ulocal[i] + extForce[i]*tau_rhoinv;
            auto usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];

            // cu(k) = c(k,i)*ueq(i), 27 3x3 dot products for D3Q27
            // not enough work to make sse/avx implementation efficient
            for (auto k=1; k<kdim; ++k){ // cu[0] = 0.0
                auto ck = &c[k*3];
                cu[k] = ck[0]*ueq[0] + ck[1]*ueq[1] + ck[2]*ueq[2];
            }
            
            // Collision - compute nprime
            auto _rholocal = _mm256_set1_ps(rholocal);
            auto nlocal = n->get(zl,yl,xl,0);
            auto _kusq = _mm256_set1_ps(-1.5f*usq); // -1.5*usq
            for (auto k=0; k<kdimAvx/8; ++k){
                auto k8 = k*8;
                auto _w = _mm256_loadu_ps(&w[k8]);
                auto _cu = _mm256_load_ps(&cu[k8]);
                auto _nk = _mm256_loadu_ps(&nlocal[k8]); // harmless access overflow
                auto _cusq = _mm256_mul_ps(_cu, _cu);
                // neq = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                auto _neq = _mm256_fmadd_ps(_three, _cu, _one); // 3*cu(k) + 1
                _neq = _mm256_fmadd_ps(_fourPointFive, _cusq, _neq); // 4.5*cusq(k) + 3*cu(k) + 1
                _neq = _mm256_add_ps(_neq, _kusq); // -1.5*usq + 4.5*cusq(k) + 3*cu(k) + 1
                _neq = _mm256_mul_ps(_neq, _rholocal);
                _neq = _mm256_mul_ps(_neq, _w);
                // nprime(k) = (1.0-_omega)*n(zl,yl,xl,k) + _omega*neq(k);
                auto _nprime = _mm256_fmadd_ps(_oneMinusOmega, _nk, _mm256_mul_ps(__omega, _neq));
                _mm256_store_ps(&nprime[k8], _nprime);
            }

            // Streaming
            for (auto k=0; k<kdim; ++k){ // NOTE: 0<=k<=kdim (NOT kdimAvx)
                auto ck = &c[k*3];
                ntmp->at(zl+ck[2],yl+ck[1],xl+ck[0],k) = nprime[k];
            }
            
        }
    }

    _mm_free(nprime);
    _mm_free(cu);
}

void BGK::_serialCollideAndStream(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    for (auto zl=1; zl<zdim+1; ++zl){
        _collideAndStreamOnPlaneAvx2(zl, lattice);
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
            _collideAndStreamOnPlaneAvx2(zl, lattice);
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
                auto nlocal = n->get(zl, yl, xl, 0);
                for (auto k=0; k<kdim; ++k){
                    rholocal += nlocal[k]; // \sum_k n[k]
                    auto ck = &c[k*3];
                    for (auto i=0; i<3; ++i)
                        ulocal[i] += nlocal[k]*ck[i]; // \sum_k n[k]*c[k][i], i=0,1,2
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
