#include "LBDynamics.hpp"

float LBDynamics::_time_collide = 0.0f;
float LBDynamics::_time_stream = 0.0f;

LBDynamics::LBDynamics(){}

LBDynamics::~LBDynamics(){}

float LBDynamics::get_time_collide() const{
    return _time_collide;
}

float LBDynamics::get_time_stream() const{
    return _time_stream;
}

float LBDynamics::get_total_time() const{
    return _time_collide+_time_stream;
}
