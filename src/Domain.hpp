#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <tuple>
#include <string>

#include "LBModel.hpp"
#include "SimData.hpp"

class Domain final {
 private:
    const LBModel *_lbmodel;
    uint32_t _xdim, _ydim, _zdim;  // domain (cuboid) dimensions
    float _fluid_viscosity;
    float _fluid_density;
    std::array<float, 3> _init_flow_velocity;
    std::array<float, 3> _external_force;

 public:
    Domain() = delete;
    // This constructor can be removed once the file
    // reading capability has been implemented
    Domain(
        const LBModel *lbmodel,
        size_t xdim, size_t ydim, size_t zdim,
        float fluid_viscosity, float fluid_density,
        std::array<float, 3> init_flow_velocity,
        std::array<float, 3> external_force);
    // Domain(std::string domain_config_file);
    Domain(Domain&) = delete;
    Domain& operator=(Domain&) = delete;
    ~Domain();
    std::tuple<size_t, size_t, size_t> get_dimensions() const;
    float get_fluid_viscosity() const;
    float get_fluid_density() const;
    float get_solid_density() const;
    const std::array<float, 3> &get_init_flow_velocity() const;
    const std::array<float, 3> &get_external_force() const;
    void initialize(SimData &simdata) const;
};

#endif
