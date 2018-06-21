#include "LBDynamics.hpp"

LBDynamics::LBDynamics(const LBModel* lbmodel, const Domain* domain){
    _collide = std::make_shared<CollisionSRT>(lbmodel, domain);
    _stream = std::make_shared<Streaming>(lbmodel, "push");
}

LBDynamics::LBDynamics(const LBModel* lbmodel, const Domain* domain, bool reference){
    _collide = std::make_shared<CollisionSRT>(lbmodel, domain, reference);
    _stream = std::make_shared<Streaming>(lbmodel, "push");
}

LBDynamics::~LBDynamics(){}

void LBDynamics::collide(SimData& simdata){
    _collide->operator()(simdata);
}

void LBDynamics::stream(SimData& simdata){
    _stream->operator()(simdata);
}

float LBDynamics::get_time_collide(){
    return _collide->get_total_time();
}

float LBDynamics::get_time_stream(){
    return _stream->get_total_time();
}

float LBDynamics::get_total_time(){
    return _collide->get_total_time() + _stream->get_total_time();
}
