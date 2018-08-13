#include "Collision.hpp"

float Collision::_time = 0.0f;

Collision::Collision() {}

Collision::~Collision() {}

float Collision::get_total_time() const {
    return _time;
}
