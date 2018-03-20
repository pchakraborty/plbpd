#include <memory>
#include <iostream>
#include "BGK.hpp"
#include "tbb/tbb.h"
#include <array>

BGK::BGK(const LBModel &lbmodel, const Domain &domain, Lattice &lattice):
    LBDynamics(), _lbmodel(lbmodel), _domain(domain), _lattice(lattice) {
    tau = 3. * _domain.getFluidViscosity() + 0.5;
    omega = 1./tau;
}

BGK::~BGK(){}

void BGK::setup(){}

void BGK::collideAndStream(){
    _parallelCollideAndStream();
}

std::array<double, 3>
BGK::_getEqlbVelocity(size_t zl, size_t yl, size_t xl,
                      size_t zdim, size_t ydim, size_t xdim){
    std::array<double, 3> ueq;
    std::array<double, 3> extForce = {0.0, 0.0, 0.0};  // TODO: Get extForce from Domain
    auto ndx3d = xdim*(ydim*zl+yl)+xl;
    for (auto i=0; i<3; ++i){
        auto ndx4d = xdim*(ydim*(zdim*i+zl)+yl)+xl;
        ueq[i] = _lattice.u[ndx4d] + extForce[i]*tau/_lattice.rho[ndx3d];
    }
    return ueq;
}

// This is the workhorse
void BGK::_collideAndStreamOnPlane(size_t zl){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDomainDimensions();

    auto numVelocityVectors = _lbmodel.getNumVelocityVectors();
    auto c = _lbmodel.getLatticeVelocities();
    auto w = _lbmodel.getDirectionalWeights();

    std::array<double, 3> extForce = {0.0, 0.0, 0.0};  // TODO: Get extForce from Domain
    for (auto yl=1; yl<ydim+1; ++yl){
        for (auto xl=1; xl<xdim+1; ++xl){
            auto ueq = _getEqlbVelocity(zl, yl, zl, zdim, ydim, xdim);
            auto usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
            auto rholoc = _lattice.rho[xdim*(ydim*zl+yl)+xl];
            for (auto k=0; k<numVelocityVectors; ++k){
                auto k3 = k*3;
                std::array<double, 3> ck = {c[k3], c[k3+1], c[k3+2]};
                auto cu = ck[0]*ueq[0] + ck[1]*ueq[1] + ck[2]*ueq[2];
                auto neq = w[k]*rholoc*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                auto old_ndx = xdim*(ydim*(zdim*k+zl)+yl)+xl;
                auto new_ndx = xdim*(ydim*(zdim*k+(zl+ck[2]))+(yl+ck[1]))+(xl+ck[0]);
                _lattice.ntmp[new_ndx] = (1.0-omega)*_lattice.n[old_ndx] + omega*neq;
            }
        }
    }
}

void BGK::_serialCollideAndStream(){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDomainDimensions();
    for (auto zl=1; zl<zdim+1; ++zl){
        _collideAndStreamOnPlane(zl);
    }
    // swap n and ntmp
    _lattice.n.swap(_lattice.ntmp);
}

void BGK::_parallelCollideAndStream(){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDomainDimensions();
    tbb::parallel_for(size_t(1), zdim+1, [=] (size_t zl){
        _collideAndStreamOnPlane(zl);
     });
    // swap n and ntmp
    _lattice.n.swap(_lattice.ntmp);
}

void BGK::calcMoments(){}

double BGK::getAvgFluidDensity(){
    return 0.0;
}

void BGK::getEqlbDist(const double rholoc, const double &uloc, double &nloc){}

void BGK::_printInfoForDebugging(){
    
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain.getDomainDimensions();
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