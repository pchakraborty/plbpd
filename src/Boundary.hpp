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

    const Domain *_domain;
    const LBModel *_lbmodel;
    size_t _xdim, _ydim, _zdim, _kdim;

    float _solidDensity;

    BoundaryType _type; // types are periodic/noslip
    BoundaryVelocity _velocity;

    bool _boundaryTypeIsPrescribed(const std::string direction) const;
    bool _boundaryVelocityIsPrescribed(const std::string direction) const;
    void _applyPeriodicityEastWest(Lattice &lattice) const;
    void _applyPeriodicityNorthSouth(Lattice &lattice) const;
    void _applyNoslipUp(Lattice &lattice) const;
    void _applyNoslipDown(Lattice &lattice) const;
    void _applyNoslipNorth(Lattice &lattice) const;
    void _applyNoslipSouth(Lattice &lattice) const;
    const std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>
    _getBoundaryExtent(const std::string direction) const;
    void _applyVelocityToBoundary(const std::string direction, Lattice &lattice) const;
    void _applyDensityToBoundary(const std::string direction, Lattice &lattice) const;

public:

    Boundary(
        const Domain *domain,
        const LBModel *lbmodel,
        const float solidDensity,
        BoundaryType type,
        BoundaryVelocity velocity
    );
    ~Boundary();
    void reset(Lattice &lattice) const;
    const BoundaryType getBoundaryType() const;
    const BoundaryVelocity getBoundaryVelocity() const;
    void applyVelocity(Lattice &lattice) const;
    void applyDensity(Lattice &lattice) const;
    void applyPeriodicity(Lattice &lattice) const;
    void applyNoslip(Lattice &lattice) const;

};

#endif
