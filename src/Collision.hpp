#ifndef COLLISION_HPP
#define COLLISION_HPP

#include <string>
#include "LBModel.hpp"
#include "Lattice.hpp"
#include <chrono>

class Collision{

protected:

    static float _time;

public:

    Collision();
    virtual ~Collision();
    Collision(Collision&) = delete;
    Collision& operator=(Collision&) = delete;
    virtual void operator()(Lattice &lattice) const = 0;
    float get_total_time() const;
};

#endif
