#ifndef LBDYNAMICS_HPP
#define LBDYNAMICS_HPP

#include <string>
#include "LBModel.hpp"
#include "Lattice.hpp"
#include <chrono>

class LBDynamics{

protected:

    static float _time_collide;
    static float _time_stream;

public:

    LBDynamics();
    virtual ~LBDynamics();
    virtual void collide(Lattice &lattice) const = 0;
    virtual void stream(Lattice &lattice) const = 0;
    float get_time_collide() const;
    float get_time_stream() const;
    float get_total_time() const;
};

#endif
