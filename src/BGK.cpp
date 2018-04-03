#include <memory>
#include <iostream>
#include "BGK.hpp"
#include "tbb/tbb.h"
#include <array>

BGK::BGK(const LBModel &lbmodel, const Domain &domain):
    LBDynamics(), _lbmodel(lbmodel), _domain(domain) {
    tau = 3. * _domain.getFluidViscosity() + 0.5;
    omega = 1./tau;
}

BGK::~BGK(){}

void BGK::initialize(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDimensions();
    auto kdim = _lbmodel.getNumVelocityVectors();
    for (auto xl=1; xl<xdim+1; ++xl){
        for (auto yl=1; yl<ydim+1; ++yl){
            for (auto zl=1; zl<zdim+1; ++zl){
                auto rholocal = _domain.getFluidDensity();
                lattice.rho[zl+(yl+xl*ydim)*zdim] = rholocal;
                auto ulocal = _domain.getInitFlowVelocity();
                for (auto i=0; i<3; ++i){
                    auto u_ndx = i+(zl+(yl+xl*ydim)*zdim)*3;
                    lattice.u[u_ndx] = ulocal[i];
                }
                std::vector<float> nlocal(kdim);
                _getEqlbDist(rholocal, ulocal, nlocal);
                for (auto k=0; k<kdim; ++k){
                    auto n_ndx = k+(zl+(yl+xl*ydim)*zdim)*kdim;
                    lattice.n[n_ndx] = nlocal[k];
                }
            }
        }
    }
}

void BGK::_getEqlbDist(const float rholocal, const std::array<float, 3> &ulocal, std::vector<float> &nlocal){
    auto usq = ulocal[0]*ulocal[0] + ulocal[1]*ulocal[1] + ulocal[2]*ulocal[2];
    auto kdim = _lbmodel.getNumVelocityVectors();
    auto w = _lbmodel.getDirectionalWeights();
    auto c = _lbmodel.getLatticeVelocities();
    for (auto k=0; k<kdim; ++k){
        auto k3 = k*3;
        std::array<int, 3> ck = {c[k3], c[k3+1], c[k3+2]};
        auto cu = ck[0]*ulocal[0] + ck[1]*ulocal[1] + ck[2]*ulocal[2];
        nlocal[k] = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    }
}

void BGK::collideAndStream(Lattice &lattice){
    _parallelCollideAndStream(lattice);
}

// This is the workhorse
void BGK::_collideAndStreamOnPlane(size_t xl, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDimensions();

    // Local variables help the compiler optimize better
    auto kdim = _lbmodel.getNumVelocityVectors();
    auto c = _lbmodel.getLatticeVelocities();
    auto w = _lbmodel.getDirectionalWeights();
    const float *rho =lattice.rho.data();
    float *n = lattice.n.data();
    float *u = lattice.u.data();
    float *ntmp = lattice.ntmp.data();

    std::array<float, 3> extForce = {0.0, 0.0, 0.0}; // TODO: Get extForce from Domain
    for (auto yl=1; yl<ydim+1; ++yl){
        for (auto zl=1; zl<zdim+1; ++zl){
            auto ndx = zl+(yl+xl*ydim)*zdim;
            auto rholocal = rho[zl+(yl+xl*ydim)*zdim];
            std::array<float, 3> ueq;
            for (auto i=0; i<3; ++i){
                auto u_ndx = i+(zl+(yl+xl*ydim)*zdim)*3;
                ueq[i] = u[u_ndx] + extForce[i]*tau/rho[ndx];
            }
            auto usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
            for (auto k=0; k<kdim; ++k){
                auto k3 = k*3;
                std::array<int, 3> ck = {c[k3], c[k3+1], c[k3+2]};
                auto cu = ck[0]*ueq[0] + ck[1]*ueq[1] + ck[2]*ueq[2];
                auto neq = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                auto n_ndx = k+(zl+(yl+xl*ydim)*zdim)*kdim;
                auto ntmp_ndx = k+((zl+ck[2])+((yl+ck[1])+(xl+ck[0])*ydim)*zdim)*kdim;
                ntmp[ntmp_ndx] = (1.0-omega)*n[n_ndx] + omega*neq;
            }
        }
    }
}

void BGK::_serialCollideAndStream(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDimensions();
    for (auto xl=1; xl<xdim+1; ++xl){
        _collideAndStreamOnPlane(xl, lattice);
    }
    // swap n and ntmp
    lattice.n.swap(lattice.ntmp);
}

void BGK::_parallelCollideAndStream(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDimensions();
    tbb::parallel_for(size_t(1), xdim+1, [this, &lattice] (size_t xl){
            _collideAndStreamOnPlane(xl, lattice);
     });
    // swap n and ntmp
    lattice.n.swap(lattice.ntmp);
}

float BGK::getAvgFluidDensity(){
    return 0.0;
}

void BGK::_printInfoForDebugging(){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDimensions();
    auto numVelocityVectors = _lbmodel.getNumVelocityVectors();
    auto c = _lbmodel.getLatticeVelocities();
    auto w = _lbmodel.getDirectionalWeights();
    
    std::cout<<"xdim: "<<xdim<<", ydim:  "<<ydim<<", zdim: "<<zdim<<std::endl;
    std::cout<<"numVelocityVectors: "<<numVelocityVectors<<std::endl;
    std::cout<<"Lattice velocities: ";
    for(auto ic=begin(c); ic!=end(c); ++ic)
        std::cout<<*ic<<" ";
    std::cout<<std::endl;
    std::cout<<"Directional weights: ";
    for(auto iw=begin(w); iw!=end(w); ++iw)
        std::cout<<*iw<<" ";
    std::cout<<std::endl;
}
