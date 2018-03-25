#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <string>
#include <unordered_map>
#include "Lattice.hpp"
#include "Domain.hpp"

using BoundaryType = std::unordered_map<std::string, std::string>;
using BoundaryVelocity = std::unordered_map<std::string, std::array<float, 3> >;

class Boundary final{

private:

    BoundaryType _type;
    BoundaryVelocity _velocity;

    void _applyEastBoundary(const Domain &domain, Lattice &lattice);
    void _applyWestBoundary(const Domain &domain, Lattice &lattice);
    void _applyNorthBoundary(const Domain &domain, Lattice &lattice);
    void _applySouthBoundary(const Domain &domain, Lattice &lattice);
    void _applyUpBoundary(const Domain &domain, Lattice &lattice);
    void _applyDownBoundary(const Domain &domain, Lattice &lattice);
    
public:

    Boundary(BoundaryType type, BoundaryVelocity velocity);
    ~Boundary();
    BoundaryType getBoundaryType() const;
    BoundaryVelocity getBoundaryVelocity() const;
    void apply(const Domain &domain, Lattice &lattice);

};

#endif
