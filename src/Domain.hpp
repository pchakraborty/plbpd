#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <string>
#include <tuple>
#include "Lattice.hpp"

class Domain final{
    
private:
    
    uint32_t _xdim, _ydim, _zdim; // domain (cuboid) dimensions
    float _fluidViscosity;
    float _fluidDensity;
    std::array<float, 3> _initFlowVelocity;

public:

    // This constructor can be removed once the file
    // reading capability has been implemented
    Domain(size_t xdim, size_t ydim, size_t zdim, float fluidViscosity, float fluidDensity, std::array<float, 3> initFlowVelocity);
    Domain(std::string domainConfigFile);
    ~Domain();
    std::tuple<size_t, size_t, size_t> getDimensions() const;
    float getFluidViscosity() const;
    float getFluidDensity() const;
    float getSolidDensity() const;
    std::array<float, 3> getInitFlowVelocity() const;
    void initialize(Lattice &lattice) const;

};

#endif
