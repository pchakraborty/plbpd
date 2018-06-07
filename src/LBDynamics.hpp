#ifndef LBDYNAMICS_HPP
#define LBDYNAMICS_HPP

#include <string>
#include "LBModel.hpp"
#include "Lattice.hpp"
#include <chrono>

class LBDynamics{

protected:

    static float _timeTakenByCollideAndStream;

public:

    LBDynamics();
    virtual ~LBDynamics();
    // virtual void collide_and_stream(Lattice &lattice) const = 0;
    virtual void collideAndStream(Lattice &lattice) const = 0; // core LB method
    float getTimeTakenByCollideAndStream() const;
    float getTotalTimeTaken() const;

};

#endif
