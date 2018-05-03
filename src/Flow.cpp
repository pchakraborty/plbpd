#include "Flow.hpp"
#include <stdexcept>

Flow::Flow(std::string flowType){
    if (flowType=="Couette2D")
        initFlowCouette2D();
    else if (flowType=="Couette3D")
        initFlowCouette3D();
    else if (flowType=="Poiseuille2D")
        initFlowPoiseuille2D();
    else if (flowType=="Poiseuille3D")
        initFlowPoiseuille3D();
    else
        throw std::invalid_argument("Flow type ["+flowType+"] was not recognized");
}

Flow::~Flow(){}

void Flow::initFlowCouette2D(){
    _lbmodel = std::make_unique<LBModel>("D2Q9");
    std::array<float, 3> initFlowVelocity = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> externalForce = {0.0f, 0.0f, 0.0f};
    _domain = std::make_unique<Domain>
        (100, 1, 25, 0.5f, 1.0f, initFlowVelocity, externalForce);
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
        (_domain.get(), _lbmodel.get(), 10.0f, type, velocity);
    _numTimeSteps = 20000;
}

void Flow::initFlowCouette3D(){
    _lbmodel = std::make_unique<LBModel>("D3Q19");
    std::array<float, 3> initFlowVelocity = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> externalForce = {0.0f, 0.0f, 0.0f};
    _domain = std::make_unique<Domain>
        (100, 25, 25, 0.5f, 1.0f, initFlowVelocity, externalForce);
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
        (_domain.get(), _lbmodel.get(), 10.0f, type, velocity);
    _numTimeSteps = 200;
}

void Flow::initFlowPoiseuille2D(){
    _lbmodel = std::make_unique<LBModel>("D2Q9");
    std::array<float, 3> initFlowVelocity = {0.0, 0.0, 0.0};
    std::array<float, 3> externalForce = {1.0e-5, 0.0, 0.0};
    _domain = std::make_unique<Domain>
        (256, 1, 48, 0.1f, 1.0f, initFlowVelocity, externalForce);
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
        (_domain.get(), _lbmodel.get(), 10.0f, type, velocity);
    _numTimeSteps = 5000;
}

void Flow::initFlowPoiseuille3D(){
    _lbmodel = std::make_unique<LBModel>("D3Q19");
    std::array<float, 3> initFlowVelocity = {0.0, 0.0, 0.0};
    std::array<float, 3> externalForce = {1.0e-5, 0.0, 0.0};
    _domain = std::make_unique<Domain>
        (256, 256, 48, 0.1f, 1.0f, initFlowVelocity, externalForce);
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
        (_domain.get(), _lbmodel.get(), 10.0f, type, velocity);
    _numTimeSteps = 10;
}

const Domain *Flow::getFlowDomain() const{
    return _domain.get();
}

const Boundary *Flow::getFlowBoundary() const{
    return _boundary.get();
}

const LBModel *Flow::getLBModel() const{
    return _lbmodel.get();
}

const size_t Flow::getNumTimeSteps() const{
    return _numTimeSteps;
}
