#include "LBDynamics.hpp"

float LBDynamics::_timeTakenByCollideAndStream = 0.0f;
float LBDynamics::_timeTakenByCalcMoments = 0.0f;

LBDynamics::LBDynamics(){}

LBDynamics::~LBDynamics(){}

float LBDynamics::getTimeTakenByCollideAndStream() const{
    return _timeTakenByCollideAndStream;
}

float LBDynamics::getTimeTakenByCalcMoments() const{
    return _timeTakenByCalcMoments;
}

float LBDynamics::getTotalTimeTaken() const{
    return _timeTakenByCollideAndStream + _timeTakenByCalcMoments;
}
