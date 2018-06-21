#include "Domain.hpp"
#include "EqlbDist.hpp"
#include <stdexcept>
#include <vector>

// Domain::Domain(std::string domain_config_file){
//     // TODO: Read domain details from domain_config_file
//     throw std::logic_error("Reading domain from config file has not been implemented");
// }

Domain::Domain(
    const LBModel *lbmodel,
    size_t xdim, size_t ydim, size_t zdim,
    float fluid_viscosity, float fluid_density,
    std::array<float, 3> init_flow_velocity,
    std::array<float, 3> external_force
               ):
    _lbmodel(lbmodel),
    _xdim(xdim), _ydim(ydim), _zdim(zdim),
    _fluid_viscosity(fluid_viscosity), _fluid_density(fluid_density),
    _init_flow_velocity(init_flow_velocity),
    _external_force(external_force)
{}
               

Domain::~Domain(){}

std::tuple<size_t, size_t, size_t> Domain::get_dimensions() const{
    return std::tie(_xdim, _ydim, _zdim);
}

float Domain::get_fluid_viscosity() const{
    return _fluid_viscosity;
}

float Domain::get_fluid_density() const{
    return _fluid_density;
}

const std::array<float, 3> &Domain::get_init_flow_velocity() const{
    return _init_flow_velocity;
}

const std::array<float, 3> &Domain::get_external_force() const{
    return _external_force;
}

void Domain::initialize(SimData &simdata) const{
    // Set rho, u, n at domain+buffer nodes, assuming that all
    // nodes are interior fluid nodes. Next, call Boundary::initialize()
    auto kdim = _lbmodel->get_num_directions();
    std::vector<float> nlocal(kdim, 0.0f);
    auto eqlbdist = EqlbDist();
    for (auto zl=0; zl<_zdim+2; ++zl){
        for (auto yl=0; yl<_ydim+2; ++yl){
            for (auto xl=0; xl<_xdim+2; ++xl){
                simdata.rho->at(zl,yl,xl) = _fluid_density;
                for (auto i=0; i<3; i++)
                    simdata.u->at(zl,yl,xl,i) = _init_flow_velocity[i];
                eqlbdist(_lbmodel, _fluid_density, _init_flow_velocity, nlocal);
                for (auto k=0; k<kdim; ++k)
                    simdata.n->at(zl,yl,xl,k) = nlocal[k];
            }
        }
    }
}
