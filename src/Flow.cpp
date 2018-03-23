#include "Flow.hpp"
#include "Domain.hpp"
#include <stdexcept>

Flow::Flow(std::string flowType){
    if (flowType=="Couette")
        initCouetteFlow();
    else
        throw std::logic_error("Flow type ["+flowType+"] was not recognized");
}

Flow::~Flow(){}

void Flow::initCouetteFlow(){
    _domain = Domain(100, 100, 100, 0.5, 1.0, 10.0, {0.0, 0.0, 0.0});
    BoundaryType type = {
        {"east", "periodic"},
        {"west", "periodic"},
        {"north", "periodic"},
        {"south", "periodic"},
        {"up", "noslip"},
        {"down", "noslip"}
    };
    BoundaryVelocity velocity = {
        {"up", {1.0, 1.0, 1.0}},
        {"down", {0.0, 0.0, 0.0}}
    };
    _boundary = std::move(std::make_unique<Boundary>(type, velocity));
}

const Domain &Flow::getFlowDomain() const{
    return _domain;
}

const Boundary *Flow::getFlowBoundary() const{
    return _boundary.get();
}
