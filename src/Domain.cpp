#include "Domain.hpp"
#include "EqlbDist.hpp"
#include <stdexcept>
#include <vector>

// Domain::Domain(std::string domainConfigFile){
//     // TODO: Read domain details from domainConfigFile
//     throw std::logic_error("Reading domain from config file has not been implemented");
// }

Domain::Domain(
    const LBModel *lbmodel,
    size_t xdim, size_t ydim, size_t zdim,
    float fluidViscosity, float fluidDensity,
    std::array<float, 3> initFlowVelocity,
    std::array<float, 3> externalForce
               ):
    _lbmodel(lbmodel),
    _xdim(xdim), _ydim(ydim), _zdim(zdim),
    _fluidViscosity(fluidViscosity), _fluidDensity(fluidDensity),
    _initFlowVelocity(initFlowVelocity),
    _externalForce(externalForce)
{}
               

Domain::~Domain(){}

std::tuple<size_t, size_t, size_t> Domain::getDimensions() const{
    return std::tie(_xdim, _ydim, _zdim);
}

float Domain::getFluidViscosity() const{
    return _fluidViscosity;
}

float Domain::getFluidDensity() const{
    return _fluidDensity;
}

const std::array<float, 3> &Domain::getInitFlowVelocity() const{
    return _initFlowVelocity;
}

const std::array<float, 3> &Domain::getExternalForce() const{
    return _externalForce;
}

void Domain::initialize(Lattice &lattice) const{
    // Set rho, u, n at domain+buffer nodes, assuming that all
    // nodes are interior fluid nodes. Next, call Boundary::initialize()
    auto kdim = _lbmodel->getNumberOfDirections();
    std::vector<float> nlocal(kdim, 0.0f);
    auto eqlbdist = EqlbDist();
    for (auto zl=0; zl<_zdim+2; ++zl){
        for (auto yl=0; yl<_ydim+2; ++yl){
            for (auto xl=0; xl<_xdim+2; ++xl){
                lattice.rho->at(zl,yl,xl) = _fluidDensity;
                for (auto i=0; i<3; i++)
                    lattice.u->at(zl,yl,xl,i) = _initFlowVelocity[i];
                eqlbdist(_lbmodel, _fluidDensity, _initFlowVelocity, nlocal);
                for (auto k=0; k<kdim; ++k)
                    lattice.n->at(zl,yl,xl,k) = nlocal[k];
            }
        }
    }
}
