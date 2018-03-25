#include "Domain.hpp"

Domain::Domain(){}

Domain::Domain(std::string domainConfigFile){
    // TODO: Read domain details from domainConfigFile
}

Domain::Domain(size_t xdim, size_t ydim, size_t zdim,
               float fluidViscosity, float fluidDensity, float solidDensity,
               std::array<float, 3> initFlowVelocity):
    _xdim(xdim), _ydim(ydim), _zdim(zdim),
    _fluidViscosity(fluidViscosity), _fluidDensity(fluidDensity),
    _solidDensity(solidDensity), _initFlowVelocity(initFlowVelocity)
{}
               

Domain::~Domain(){}

std::tuple<size_t, size_t, size_t> Domain::getDomainDimensions() const{
    return std::tie(_xdim, _ydim, _zdim);
}

float Domain::getFluidViscosity() const{
    return _fluidViscosity;
}

float Domain::getFluidDensity() const{
    return _fluidDensity;
}

float Domain::getSolidDensity() const{
    return _solidDensity;
}

std::array<float, 3> Domain::getInitFlowVelocity() const{
    return _initFlowVelocity;
}
