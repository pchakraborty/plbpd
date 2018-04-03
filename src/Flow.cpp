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
    _domain = Domain(100, 100, 100, 0.5f, 1.0f, 10.0f, {0.0f, 0.0f, 0.0f});
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
    _boundary = std::move(std::make_unique<Boundary>(_domain, type, velocity));
}

const Domain &Flow::getFlowDomain() const{
    return _domain;
}

Boundary *Flow::getFlowBoundary() const{
    return _boundary.get();
}
