#include "LBDynamics.hpp"

#include <string>

LBDynamics::LBDynamics(const LBModel* lbmodel, const Domain* domain) {
    _collide = std::make_unique<CollisionSRT>(lbmodel, domain);
    _stream = std::make_unique<Streaming>(lbmodel);
}

LBDynamics::~LBDynamics() {}

void LBDynamics::collide(SimData& simdata, bool reference) const {
    _collide->operator()(simdata, reference);
}

void LBDynamics::stream(
    SimData& simdata,
    std::string stream_type,
    bool reference) const {
    _stream->operator()(simdata, stream_type, reference);
}

float LBDynamics::get_time_collide() const {
    return _collide->get_total_time();
}

float LBDynamics::get_time_stream() const {
    return _stream->get_total_time();
}

float LBDynamics::get_total_time() const {
    return _collide->get_total_time() + _stream->get_total_time();
}
