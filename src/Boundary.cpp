#include "Boundary.hpp"

#include <chrono>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include "tbb/tbb.h"

namespace chrono = std::chrono;

float Boundary::_time_noslip = 0.0f;
float Boundary::_time_periodicity = 0.0f;
float Boundary::_time_reset = 0.0f;
float Boundary::_time_velocity_reset = 0.0f;
float Boundary::_time_density_reset = 0.0f;
const std::vector<std::string> Boundary::_directions =
    {"east", "west", "north", "south", "up", "down"};

Boundary::Boundary(
    const LBModel *lbmodel, const Domain *domain,
    float solid_density, BoundaryType type, BoundaryVelocity velocity)
    : _kdim(lbmodel->get_num_directions()),
      _c(lbmodel->get_directional_velocities()),
      _w(lbmodel->get_directional_weights()),
      _reverse(lbmodel->get_reverse()),
      _cs2inv(1.0f/lbmodel->get_speed_of_sound_squared()),
      _solid_density(solid_density),
      _type(type),
      _velocity(velocity) {
    std::tie(_xdim, _ydim, _zdim) = domain->get_dimensions();
}

Boundary::~Boundary() {}

const BoundaryType Boundary::get_boundary_type() const {
    return _type;
}

const BoundaryVelocity Boundary::get_boundary_velocity() const {
    return _velocity;
}

bool Boundary::_velocity_is_prescribed(std::string direction) const {
    if (_velocity.find(direction) != _velocity.end())
        return true;
    else
        return false;
}

bool Boundary::_type_is_prescribed(std::string direction) const {
    if (_type.find(direction) != _type.end())
        return true;
    else
        return false;
}

void Boundary::reset(SimData &simdata) const {
    auto start = chrono::system_clock::now();

    reset_velocity(simdata);
    reset_density(simdata);

    chrono::duration<float> elapsed = chrono::system_clock::now()-start;
    Boundary::_time_reset += elapsed.count();
}

void Boundary::reset_velocity(SimData &simdata) const {
    auto start = chrono::system_clock::now();

    for (const std::string& dirxn : Boundary::_directions)
        if (_velocity_is_prescribed(dirxn))
            _set_velocity(dirxn, simdata.u);

    chrono::duration<float> elapsed = chrono::system_clock::now()-start;
    Boundary::_time_velocity_reset += elapsed.count();
}

void Boundary::reset_density(SimData &simdata) const {
    auto start = chrono::system_clock::now();

    for (const std::string& dirxn : Boundary::_directions)
        if (_type_is_prescribed(dirxn))
            if (_type.at(dirxn) == "noslip")
                _set_density(dirxn, simdata.rho);

    chrono::duration<float> elapsed = chrono::system_clock::now()-start;
    Boundary::_time_density_reset += elapsed.count();
}

void Boundary::apply_periodicity(SimData &simdata) const {
    auto start = chrono::system_clock::now();

    _apply_periodicity_east_west(simdata);
    _apply_periodicity_north_south(simdata);

    chrono::duration<float> elapsed = chrono::system_clock::now()-start;
    Boundary::_time_periodicity += elapsed.count();
}

void Boundary::apply_noslip(SimData &simdata) const {
    auto start = chrono::system_clock::now();

    for (const std::string& dirxn : Boundary::_directions)
        if (_type_is_prescribed(dirxn))
            if (_type.at(dirxn) == "noslip")
                _apply_noslip(dirxn, simdata);

    chrono::duration<float> elapsed = chrono::system_clock::now()-start;
    Boundary::_time_noslip += elapsed.count();
}

const std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>
Boundary::_get_boundary_extent(std::string direction) const {
    if (direction == "up") {
        const uint32_t zmin = _zdim, zmax = _zdim + 1;
        const uint32_t ymin = 1, ymax = _ydim + 1;
        const uint32_t xmin = 1, xmax = _xdim + 1;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction == "down") {
        const uint32_t zmin = 1, zmax = 1 + 1;
        const uint32_t ymin = 1, ymax = _ydim + 1;
        const uint32_t xmin = 1, xmax = _xdim + 1;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction == "north") {
        const uint32_t zmin = 1, zmax = _zdim + 1;
        const uint32_t ymin = _ydim, ymax = _ydim + 1;
        const uint32_t xmin = 1, xmax = _xdim + 1;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction == "south") {
        const uint32_t zmin = 1, zmax = _zdim + 1;
        const uint32_t ymin = 1, ymax = 1 + 1;
        const uint32_t xmin = 1, xmax = _xdim + 1;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction == "east") {
        const uint32_t zmin = 1, zmax = _zdim + 1;
        const uint32_t ymin = 1, ymax = _ydim + 1;
        const uint32_t xmin = _xdim, xmax = _xdim + 1;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction == "west") {
        const uint32_t zmin = 1, zmax = _zdim + 1;
        const uint32_t ymin = 1, ymax = _ydim + 1;
        const uint32_t xmin = 1, xmax = 1 + 1;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else {
        throw std::logic_error("invalid direction " + direction);
    }
}

void Boundary::_set_velocity(std::string direction, fVectorField *u) const {
    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin,
             ymax, zmin, zmax) = _get_boundary_extent(direction);

    tbb::parallel_for
    (tbb::blocked_range3d<uint32_t> (zmin, zmax, ymin, ymax, xmin, xmax),
     [this, u, &direction]
     (const tbb::blocked_range3d<uint32_t> &r) {
        // lambda body - start
        for (auto zl = r.pages().begin(); zl < r.pages().end(); ++zl)
            for (auto yl = r.rows().begin(); yl < r.rows().end(); ++yl)
                for (auto xl = r.cols().begin(); xl < r.cols().end(); ++xl) {
                    auto ulocal = u->get(zl, yl, xl, 0);
                    for (auto i = 0; i < 3; ++i)
                        ulocal[i] = _velocity.at(direction)[i];
                }
        // lambda body - end
    });
}

void Boundary::_set_density(std::string direction, fScalarField *rho) const {
    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin,
             ymax, zmin, zmax) = _get_boundary_extent(direction);

    tbb::parallel_for
    (tbb::blocked_range3d<uint32_t> (zmin, zmax, ymin, ymax, xmin, xmax),
     [this, rho]
     (const tbb::blocked_range3d<uint32_t> &r) {
        // lambda body - start
        for (auto zl = r.pages().begin(); zl < r.pages().end(); ++zl)
            for (auto yl = r.rows().begin(); yl < r.rows().end(); ++yl)
                for (auto xl = r.cols().begin(); xl < r.cols().end(); ++xl)
                    rho->at(zl, yl, xl) = _solid_density;
        // lambda body - end
    });
}

void Boundary::_apply_noslip(std::string direction, SimData &simdata) const {
    auto n = simdata.n;
    auto rho = simdata.rho;

    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin,
             ymax, zmin, zmax) = _get_boundary_extent(direction);

    const auto ub = _velocity.at(direction);

    tbb::parallel_for
    (tbb::blocked_range3d<uint32_t> (zmin, zmax, ymin, ymax, xmin, xmax),
     [this, &ub, n, rho]
     (const tbb::blocked_range3d<uint32_t> &r) {
        // lambda body - start
        for (auto zl = r.pages().begin(); zl < r.pages().end(); ++zl)
            for (auto yl = r.rows().begin(); yl < r.rows().end(); ++yl)
                for (auto xl = r.cols().begin(); xl < r.cols().end(); ++xl)
                    // kp = 0 => current node
                    for (auto kp = 1; kp < _kdim; ++kp) {
                        size_t nz, ny, nx;  // kp-th neighbor of (zl, yl, xl)
                        std::tie(nz, ny, nx) =
                            n->get_neighbor(zl, yl, xl, _c[kp]);
                        auto rhonbr = rho->at(nz, ny, nx);
                        auto k = _reverse[kp];
                        auto ck = _c[k].data();
                        auto cu = ck[0]*ub[0] + ck[1]*ub[1] + ck[2]*ub[2];
                        n->at(nz, ny, nx, kp) =
                            n->at(zl, yl, xl, k) - 2.0f*_w[k]*rhonbr*cu*_cs2inv;
                    }
        // lambda body -end
    });
}

bool Boundary::_is_east_west_periodic() const {
    bool periodic = false;
    if (_type_is_prescribed("east") && _type_is_prescribed("west")) {
        if ( (_type.at("east") == "periodic") &&
             (_type.at("west") == "periodic") ) {
            periodic = true;
        }
    }
    return periodic;
}

void Boundary::_apply_periodicity_east_west(SimData &simdata) const {
    if (_is_east_west_periodic()) {
        auto n = simdata.n;
        for (auto zl = 1; zl < _zdim+1; ++zl) {
            for (auto yl = 1; yl < _ydim+1; ++yl) {
                for (auto k = 1; k < _kdim; ++k) {  // k=0 => rest particle
                    // east->west (x-positive components)
                    auto ck0 = _c[k][0];
                    if (ck0 > 0)
                        n->at(zl, yl, 1, k) = n->at(zl, yl, _xdim+1, k);
                    // west->east (x-negative components)
                    if (ck0 < 0)
                        n->at(zl, yl, _xdim, k) = n->at(zl, yl, 0, k);
                }
            }
        }
    }
}

bool Boundary::_is_north_south_periodic() const {
    bool periodic = false;
    if (_type_is_prescribed("north") && _type_is_prescribed("south")) {
        if ( (_type.at("north") == "periodic") &&
             (_type.at("south") == "periodic") ) {
            periodic = true;
        }
    }
    return periodic;
}

void Boundary::_apply_periodicity_north_south(SimData &simdata) const {
    if (_is_north_south_periodic()) {
        auto n = simdata.n;
        for (auto zl = 1; zl < _zdim+1; ++zl) {
            for (auto xl = 1; xl < _xdim+1; ++xl) {
                for (auto k = 1; k < _kdim; ++k) {  // k=0 => rest particle
                    auto ck1 = _c[k][1];
                    // north->south (y-positive components)
                    if (ck1 > 0)
                        n->at(zl, 1, xl, k) = n->at(zl, _ydim+1, xl, k);
                    // south->north (y-negative components)
                    if (ck1 < 0)
                        n->at(zl, _ydim, xl, k) = n->at(zl, 0, xl, k);
                }
            }
        }
    }
}

float Boundary::get_time_noslip() const {
    return _time_noslip;
}

float Boundary::get_time_periodicity() const {
    return _time_periodicity;
}

float Boundary::get_time_reset() const {
    return _time_reset;
}

float Boundary::get_time_velocity_reset() const {
    return _time_velocity_reset;
}

float Boundary::get_time_density_reset() const {
    return _time_density_reset;
}

float Boundary::get_total_time() const {
    return _time_noslip + _time_periodicity + _time_reset;
}
