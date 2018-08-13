#ifndef COLLISION_HPP
#define COLLISION_HPP

#include <chrono>
#include <string>

#include "LBModel.hpp"
#include "SimData.hpp"

class Collision {
 protected:
    static float _time;

 public:
    Collision();
    Collision(Collision&) = delete;
    Collision& operator=(Collision&) = delete;
    virtual ~Collision();
    virtual void operator()(SimData &simdata, bool reference) const = 0;
    float get_total_time() const;
};

#endif
