#include "Flow.hpp"
#include <stdexcept>

Flow::Flow(std::string flowType){
    if (flowType=="Couette")
        initCouetteFlow();
    else
        throw std::logic_error("Flow type ["+flowType+"] was not recognized");
}

Flow::~Flow(){}

void Flow::initCouetteFlow(){
    std::array<float, 3> initFlowVelocity = {0.0f, 0.0f, 0.0f};
    _domain = std::make_unique<Domain>(1000, 10, 100, 0.5f, 1.0f, 10.0f, initFlowVelocity);
    BoundaryType type = {
        {"east", "periodic"},
        {"west", "periodic"},
        {"north", "periodic"},
        {"south", "periodic"},
        {"up", "noslip"},
        {"down", "noslip"}
    };
    BoundaryVelocity velocity = {
        {"up", {1.0f, 1.0f, 1.0f}},
        {"down", {0.0f, 0.0f, 0.0f}}
    };
    _boundary = std::make_unique<Boundary>(_domain.get(), type, velocity);
}

const Domain *Flow::getFlowDomain() const{
    return _domain.get();
}

const Boundary *Flow::getFlowBoundary() const{
    return _boundary.get();
}
