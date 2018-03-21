#include "Domain.hpp"

Domain::Domain(std::string domainConfigFile){
    // TODO: Read dimensions from domainConfigFile
    _xdim = 100;
    _ydim = 100;
    _zdim = 100; // 1 => 2D
    _fluidViscosity = 0.2;
}

std::tuple<size_t, size_t, size_t> Domain::getDomainDimensions() const{
    return std::tie(_xdim, _ydim, _zdim);
}

float Domain::getFluidViscosity() const{
    return _fluidViscosity;
}

Domain::~Domain(){}
