#include "Flow.hpp"
#include <stdexcept>

Flow::Flow(std::string flow_type){
    if (flow_type=="Couette2D")
        init_flow_couette_2d();
    else if (flow_type=="Couette3D")
        init_flow_couette_3d();
    else if (flow_type=="Poiseuille2D")
        init_flow_poiseuille_2d();
    else if (flow_type=="Poiseuille3D")
        init_flow_poiseuille_3d();
    else
        throw std::invalid_argument("Flow type ["+flow_type+"] was not recognized");
}

Flow::~Flow(){}

void Flow::init_flow_couette_2d(){
    _lbmodel = std::make_unique<LBModel>("D2Q9");
    std::array<float, 3> init_flow_velocity = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> ext_force = {0.0f, 0.0f, 0.0f};
    _domain = std::make_unique<Domain>
        (_lbmodel.get(), 100, 1, 25, 0.5f, 1.0f, init_flow_velocity, ext_force);
    BoundaryType type = {
        {"east", "periodic"},
        {"west", "periodic"},
        {"up", "noslip"},
        {"down", "noslip"}
    };
    BoundaryVelocity velocity = {
        {"up", {0.5f, 0.0f, 0.0f}},
        {"down", {0.0f, 0.0f, 0.0f}}
    };
    _boundary = std::make_unique<Boundary>
        (_lbmodel.get(), _domain.get(), 1.0f, type, velocity);
    _num_timesteps = 20000;
}

void Flow::init_flow_couette_3d(){
    _lbmodel = std::make_unique<LBModel>("D3Q19");
    std::array<float, 3> init_flow_velocity = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> ext_force = {0.0f, 0.0f, 0.0f};
    _domain = std::make_unique<Domain>
        (_lbmodel.get(), 100, 25, 25, 0.5f, 1.0f, init_flow_velocity, ext_force);
    BoundaryType type = {
        {"east", "periodic"},
        {"west", "periodic"},
        {"north", "periodic"},
        {"south", "periodic"},
        {"up", "noslip"},
        {"down", "noslip"}
    };
    BoundaryVelocity velocity = {
        {"up", {0.5f, 0.0f, 0.0f}},
        {"down", {0.0f, 0.0f, 0.0f}}
    };
    _boundary = std::make_unique<Boundary>
        (_lbmodel.get(), _domain.get(), 1.0f, type, velocity);
    _num_timesteps = 200;
}

void Flow::init_flow_poiseuille_2d(){
    _lbmodel = std::make_unique<LBModel>("D2Q9");
    std::array<float, 3> init_flow_velocity = {0.0, 0.0, 0.0};
    std::array<float, 3> ext_force = {1.0e-5, 0.0, 0.0};
    _domain = std::make_unique<Domain>
        (_lbmodel.get(), 256, 1, 48, 0.1f, 1.0f, init_flow_velocity, ext_force);
    BoundaryType type = {
        {"east", "periodic"},
        {"west", "periodic"},
        {"up", "noslip"},
        {"down", "noslip"}
    };
    BoundaryVelocity velocity = {
        {"up", {0.0f, 0.0f, 0.0f}},
        {"down", {0.0f, 0.0f, 0.0f}}
    };
    _boundary = std::make_unique<Boundary>
        (_lbmodel.get(), _domain.get(), 1.0f, type, velocity);
    _num_timesteps = 5000;
}

void Flow::init_flow_poiseuille_3d(){
    _lbmodel = std::make_unique<LBModel>("D3Q19");
    std::array<float, 3> init_flow_velocity = {0.0, 0.0, 0.0};
    std::array<float, 3> ext_force = {1.0e-5, 0.0, 0.0};
    _domain = std::make_unique<Domain>
        (_lbmodel.get(), 256, 256, 48, 0.1f, 1.0f, init_flow_velocity, ext_force);
    BoundaryType type = {
        {"east", "periodic"},
        {"west", "periodic"},
        {"north", "periodic"},
        {"south", "periodic"},
        {"up", "noslip"},
        {"down", "noslip"}
    };
    BoundaryVelocity velocity = {
        {"up", {0.0f, 0.0f, 0.0f}},
        {"down", {0.0f, 0.0f, 0.0f}}
    };
    _boundary = std::make_unique<Boundary>
        (_lbmodel.get(), _domain.get(), 1.0f, type, velocity);
    _num_timesteps = 10;
}

const Domain *Flow::get_domain() const{
    return _domain.get();
}

const Boundary *Flow::get_boundary() const{
    return _boundary.get();
}

const LBModel *Flow::get_lbmodel() const{
    return _lbmodel.get();
}

size_t Flow::get_num_timesteps() const{
    return _num_timesteps;
}
