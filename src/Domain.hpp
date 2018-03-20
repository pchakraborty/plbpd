#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <string>
#include <vector>
#include <tuple>

enum class NodeType {INTFLD, INFLOW, OUTFLOW, NOSLIP, FREESLIP};

class Domain final{
    
private:
    
    size_t _xdim, _ydim, _zdim; // domain (cuboid) dimensions
    double _fluidViscosity;
    
public:
    
    Domain(std::string domainConfigFile);
    ~Domain();
    std::tuple<size_t, size_t, size_t> getDomainDimensions() const;
    // getBoundaryInformation();
    double getFluidViscosity() const;
    // void applyBoundaryConditions();

};

#endif
