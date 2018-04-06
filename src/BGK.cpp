#include <memory>
#include <iostream>
#include "BGK.hpp"
#include "tbb/tbb.h"
#include <array>

BGK::BGK(const LBModel *lbmodel, const Domain *domain):
    LBDynamics(), _lbmodel(lbmodel), _domain(domain) {
    tau = 3. * _domain->getFluidViscosity() + 0.5;
    omega = 1./tau;
}

BGK::~BGK(){}

void BGK::initialize(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    auto kdim = _lbmodel->getNumVelocityVectors();
    for (auto zl=1; zl<zdim+1; ++zl){
        for (auto yl=1; yl<ydim+1; ++yl){
            auto tmp = (yl+zl*(ydim+2))*(xdim+2);
            for (auto xl=1; xl<xdim+1; ++xl){
                auto zyx = xl+tmp;
                auto rholocal = _domain->getFluidDensity();
                lattice.rho[zyx] = rholocal;
                auto ulocal = _domain->getInitFlowVelocity();
                auto zyxTimes3 = zyx*3;
                for (auto i=0; i<3; ++i)
                    lattice.u[i+zyxTimes3] = ulocal[i];
                std::vector<float> nlocal(kdim);
                _getEqlbDist(rholocal, ulocal, nlocal);
                auto zyxTimesKdim = zyx*kdim;
                for (auto k=0; k<kdim; ++k)
                    lattice.n[k+zyxTimesKdim] = nlocal[k];
            }
        }
    }
}

void BGK::_getEqlbDist(const float rholocal, const std::array<float, 3> &ulocal, std::vector<float> &nlocal){
    auto usq = ulocal[0]*ulocal[0] + ulocal[1]*ulocal[1] + ulocal[2]*ulocal[2];
    auto kdim = _lbmodel->getNumVelocityVectors();
    auto w = _lbmodel->getDirectionalWeights();
    auto c = _lbmodel->getLatticeVelocities();
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
void BGK::_collideAndStreamOnPlane(size_t zl, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();

    const auto kdim = _lbmodel->getNumVelocityVectors();
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeights();

    const float *rho = lattice.rho.data();
    const float *u = lattice.u.data();
    float *n = lattice.n.data();
    float *ntmp = lattice.ntmp.data();

    // const float *rho = lattice.rho;
    // const float *u = lattice.u;
    // float *n = lattice.n;
    // float *ntmp = lattice.ntmp;
    
    std::array<float, 3> extForce = {0.0, 0.0, 0.0}; // TODO: Get extForce from Domain
    
    for (auto yl=1; yl<ydim+1; ++yl){
        auto tmp = (yl+zl*(ydim+2))*(xdim+2);
        for (auto xl=1; xl<xdim+1; ++xl){
            auto zyx = xl+tmp;
            auto rholocal = rho[zyx];
            std::array<float, 3> ueq;
            auto zyxTimes3 = zyx*3;
            for (auto i=0; i<3; ++i)
                ueq[i] = u[i+zyxTimes3] + extForce[i]*tau/rholocal;
            auto usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
            auto zyxTimesKdim = zyx*kdim;
            for (auto k=0; k<kdim; ++k){
                auto k3 = k*3;
                std::array<int, 3> ck = {c[k3], c[k3+1], c[k3+2]};
                auto cu = ck[0]*ueq[0] + ck[1]*ueq[1] + ck[2]*ueq[2];
                auto neq = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                auto zyxktmp =  k+((xl+ck[0])+((yl+ck[1])+(zl+ck[2])*(ydim+2))*(zdim+2))*kdim;
                ntmp[zyxktmp] = (1.0-omega)*n[k+zyxTimesKdim] + omega*neq;
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
    lattice.n.swap(lattice.ntmp);
    // float *tmp = lattice.n;
    // lattice.n = lattice.ntmp;
    // lattice.ntmp = tmp;
}

void BGK::_parallelCollideAndStream(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    tbb::parallel_for(size_t(1), zdim+1, [this, &lattice] (size_t zl){
            _collideAndStreamOnPlane(zl, lattice);
     });
    // swap n and ntmp
    lattice.n.swap(lattice.ntmp);
    // float *tmp = lattice.n;
    // lattice.n = lattice.ntmp;
    // lattice.ntmp = tmp;
}

float BGK::getAvgFluidDensity(){
    return 0.0;
}

void BGK::calcMoments(Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    
    const auto kdim = _lbmodel->getNumVelocityVectors();
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeights();

    const float *n = lattice.n.data();
    float *rho =lattice.rho.data();
    float *u = lattice.u.data();

    // const float *n = lattice.n;
    // float *rho =lattice.rho;
    // float *u = lattice.u;

    for (auto zl=1; zl<zdim+1; ++zl){
        for (auto yl=1; yl<ydim+1; ++yl){
            auto tmp = (yl+zl*(ydim+2))*(xdim+2);
            for (auto xl=1; xl<xdim+1; ++xl){
                float rholocal = 0.0f;
                std::array<float, 3> ulocal = {0.0f, 0.0f, 0.0f};
                auto zyx = xl+tmp;
                auto zyxTimesKdim = zyx*kdim;
                for (auto k=0; k<kdim; ++k){
                    auto nk = n[k+zyxTimesKdim];
                    rholocal += nk;
                    auto k3 = k*3;
                    std::array<int, 3> ck = {c[k3], c[k3+1], c[k3+2]};
                    ulocal[0] += nk*ck[0];
                    ulocal[1] += nk*ck[1];
                    ulocal[2] += nk*ck[2];
                }
                rho[zyx] = rholocal;
                auto rhoinv = 1.0f/rholocal;
                auto zyxTimes3 = zyx*3;
                for (auto i=0; i<3; ++i)
                    u[i+zyxTimes3] = ulocal[i];
            }
        }
    }
}

void BGK::_printInfoForDebugging(){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    auto numVelocityVectors = _lbmodel->getNumVelocityVectors();
    auto c = _lbmodel->getLatticeVelocities();
    auto w = _lbmodel->getDirectionalWeights();
    
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
