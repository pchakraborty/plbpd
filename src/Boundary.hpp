#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <string>
#include <memory>

#include <unordered_map>
#include "Lattice.hpp"
#include "Domain.hpp"

using BoundaryType = std::unordered_map<std::string, std::string>;
using BoundaryVelocity = std::unordered_map<std::string, std::array<float, 3> >;

class Boundary final{

private:

    size_t _xdim, _ydim, _zdim;
    BoundaryType _type;
    BoundaryVelocity _velocity;

    void _applyEastBoundary(Lattice &lattice);
    void _applyWestBoundary(Lattice &lattice);
    void _applyNorthBoundary(Lattice &lattice);
    void _applySouthBoundary(Lattice &lattice);
    void _applyUpBoundary(Lattice &lattice);
    void _applyDownBoundary(Lattice &lattice);
    
public:

    Boundary(
        const Domain &domain,
        BoundaryType type,
        BoundaryVelocity velocity
    );
    ~Boundary();
    BoundaryType getBoundaryType() const;
    BoundaryVelocity getBoundaryVelocity() const;
    void apply(Lattice &lattice);

};

#endif
