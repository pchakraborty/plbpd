#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <string>
#include <memory>
#include <unordered_map>

#include "SimData.hpp"
#include "Domain.hpp"
#include "LBModel.hpp"

using BoundaryType = std::unordered_map<std::string, std::string>;
using BoundaryVelocity = std::unordered_map<std::string, std::array<float, 3> >;

class Boundary final {
 private:
    const Domain *_domain;
    const LBModel *_lbmodel;
    size_t _xdim, _ydim, _zdim, _kdim;

    float _solid_density;

    BoundaryType _type;  // types are periodic/noslip
    BoundaryVelocity _velocity;

    bool _boundary_type_is_prescribed(const std::string direction) const;
    bool _boundary_velocity_is_prescribed(const std::string direction) const;
    void _apply_periodicity_east_west(SimData &simdata) const;
    void _apply_periodicity_north_south(SimData &simdata) const;
    const std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>
    _get_boundary_extent(const std::string direction) const;
    void _apply_velocity_to_boundary(const std::string direction, SimData &simdata) const;
    void _apply_density_to_boundary(const std::string direction, SimData &simdata) const;
    void _apply_noslip_to_boundary(const std::string direction, SimData &simdata) const;

    // Timers
    static float _time_noslip;
    static float _time_periodicity;
    static float _time_reset;
    static float _time_velocity_reset;
    static float _time_density_reset;

    // List of directions (north, south etc.)
    static const std::vector<std::string> _directions;

 public:
    Boundary() = delete;
    Boundary(
        const LBModel *lbmodel,
        const Domain *domain,
        const float solid_density,
        BoundaryType type,
        BoundaryVelocity velocity);
    Boundary(Boundary&) = delete;
    Boundary& operator=(Boundary&) = delete;
    ~Boundary();

    void reset(SimData &simdata) const;
    const BoundaryType get_boundary_type() const;
    const BoundaryVelocity get_boundary_velocity() const;
    void apply_velocity(SimData &simdata) const;
    void apply_density(SimData &simdata) const;
    void apply_periodicity(SimData &simdata) const;
    void apply_noslip(SimData &simdata) const;
    // Timer access
    float get_time_noslip() const;
    float get_time_periodicity() const;
    float get_time_reset() const;
    float get_time_velocity_reset() const;
    float get_time_density_reset() const;
    float get_total_time() const;
};

#endif
