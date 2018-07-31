#include "LBDynamics.hpp"

void LBDynamics::_init(
    const LBModel* lbmodel,
    const Domain* domain,
    std::string stream_type,
    bool reference){

    _collide = std::make_shared<CollisionSRT>(lbmodel, domain, reference);
    _stream = std::make_shared<Streaming>(lbmodel, stream_type, reference);

}

LBDynamics::LBDynamics(const LBModel* lbmodel, const Domain* domain){
    _init(lbmodel, domain, "push", false);
}

LBDynamics::LBDynamics(const LBModel* lbmodel, const Domain* domain, bool reference){
    _init(lbmodel, domain, "push", reference);
}

LBDynamics::LBDynamics(const LBModel* lbmodel, const Domain* domain, std::string stream_type){
    _init(lbmodel, domain, stream_type, false);
}

LBDynamics::LBDynamics(
    const LBModel* lbmodel,
    const Domain* domain,
    std::string stream_type,
    bool reference){

    _init(lbmodel, domain, stream_type, reference);

}

LBDynamics::~LBDynamics(){}

void LBDynamics::collide(SimData& simdata) const{
    _collide->operator()(simdata);
}

void LBDynamics::stream(SimData& simdata) const{
    _stream->operator()(simdata);
}

float LBDynamics::get_time_collide() const{
    return _collide->get_total_time();
}

float LBDynamics::get_time_stream() const{
    return _stream->get_total_time();
}

float LBDynamics::get_total_time() const{
    return _collide->get_total_time() + _stream->get_total_time();
}
