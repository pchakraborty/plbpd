#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <string>
#include <memory>

#include <unordered_map>
#include "Lattice.hpp"
#include "Domain.hpp"
#include "LBModel.hpp"

using BoundaryType = std::unordered_map<std::string, std::string>;
using BoundaryVelocity = std::unordered_map<std::string, std::array<float, 3> >;

class Boundary final{

private:

    size_t _xdim, _ydim, _zdim, _kdim;
    BoundaryType _type;
    BoundaryVelocity _velocity;

    bool _boundaryTypeIsPrescribed(const std::string direction) const;
    bool _boundaryVelocityIsPrescribed(const std::string direction) const;
    void _applyEastBoundary(Lattice &lattice) const;
    void _applyWestBoundary(Lattice &lattice) const;
    void _applyNorthBoundary(Lattice &lattice) const;
    void _applySouthBoundary(Lattice &lattice) const;
    void _applyUpBoundary(Lattice &lattice) const;
    void _applyDownBoundary(Lattice &lattice) const;
    
public:

    Boundary(
        const Domain *domain,
        const LBModel *lbmodel,
        BoundaryType type,
        BoundaryVelocity velocity
    );
    ~Boundary();
    const BoundaryType getBoundaryType() const;
    const BoundaryVelocity getBoundaryVelocity() const;
    void apply(Lattice &lattice) const;

};

#endif
