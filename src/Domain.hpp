#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <string>
#include <tuple>
#include "LBModel.hpp"
#include "Lattice.hpp"

class Domain final{
    
private:
    
    const LBModel *_lbmodel;
    uint32_t _xdim, _ydim, _zdim; // domain (cuboid) dimensions
    float _fluidViscosity;
    float _fluidDensity;
    std::array<float, 3> _initFlowVelocity;
    std::array<float, 3> _externalForce;

public:

    // This constructor can be removed once the file
    // reading capability has been implemented
    Domain(
        const LBModel *lbmodel,
        size_t xdim, size_t ydim, size_t zdim,
        float fluidViscosity, float fluidDensity,
        std::array<float, 3> initFlowVelocity,
        std::array<float, 3> externalForce
    );
    //Domain(std::string domainConfigFile);
    ~Domain();
    std::tuple<size_t, size_t, size_t> getDimensions() const;
    float getFluidViscosity() const;
    float getFluidDensity() const;
    float getSolidDensity() const;
    const std::array<float, 3> &getInitFlowVelocity() const;
    const std::array<float, 3> &getExternalForce() const;
    void initialize(Lattice &lattice) const;

};

#endif
