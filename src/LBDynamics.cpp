#include "LBDynamics.hpp"

float LBDynamics::_timeTakenByCollideAndStream = 0.0f;

LBDynamics::LBDynamics(){}

LBDynamics::~LBDynamics(){}

float LBDynamics::getTimeTakenByCollideAndStream() const{
    return _timeTakenByCollideAndStream;
}

float LBDynamics::getTotalTimeTaken() const{
    return _timeTakenByCollideAndStream;
}
