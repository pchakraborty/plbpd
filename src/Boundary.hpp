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

    float _solid_density;

    BoundaryType _type; // types are periodic/noslip
    BoundaryVelocity _velocity;

    bool _boundary_type_is_prescribed(const std::string direction) const;
    bool _boundary_velocity_is_prescribed(const std::string direction) const;
    void _apply_periodicity_east_west(Lattice &lattice) const;
    void _apply_periodicity_north_south(Lattice &lattice) const;
    const std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>
    _get_boundary_extent(const std::string direction) const;
    void _apply_velocity_to_boundary(const std::string direction, Lattice &lattice) const;
    void _apply_density_to_boundary(const std::string direction, Lattice &lattice) const;
    void _apply_noslip_to_boundary(const std::string direction, Lattice &lattice) const;

    // Timers
    static float _time_noslip;
    static float _time_periodicity;
    static float _time_reset;
    
public:

    Boundary(
        const LBModel *lbmodel,
        const Domain *domain,
        const float solid_density,
        BoundaryType type,
        BoundaryVelocity velocity
    );
    ~Boundary();
    void reset(Lattice &lattice) const;
    const BoundaryType get_boundary_type() const;
    const BoundaryVelocity get_boundary_velocity() const;
    void apply_velocity(Lattice &lattice) const;
    void apply_density(Lattice &lattice) const;
    void apply_periodicity(Lattice &lattice) const;
    void apply_noslip(Lattice &lattice) const;
    // Timer access
    float get_time_noslip() const;
    float get_time_periodicity() const;
    float get_time_reset() const;
    float get_total_time() const;
    
};

#endif
