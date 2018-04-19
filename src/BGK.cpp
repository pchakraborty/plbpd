#include <memory>
#include <iostream>
#include "BGK.hpp"
#include "tbb/tbb.h"
#include <array>

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
                auto k3 = k*3;
                auto cu = c[k3+0]*ueq[0] + c[k3+1]*ueq[1] + c[k3+2]*ueq[2];
                auto neq = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                ntmp->at(zl+c[k3+2],yl+c[k3+1],xl+c[k3+0],k) =
                    (1.0f-_omega)*n->at(zl,yl,xl,k) + _omega*neq;
            }
        }
    }
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
