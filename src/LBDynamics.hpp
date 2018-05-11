#ifndef LBDYNAMICS_HPP
#define LBDYNAMICS_HPP

#include <string>
#include "LBModel.hpp"
#include "Lattice.hpp"
#include <chrono>

class LBDynamics{

protected:

    static float _timeTakenByCollideAndStream;
    static float _timeTakenByCalcMoments;

public:

    LBDynamics();
    virtual ~LBDynamics();
    virtual void collideAndStream(Lattice &lattice) const = 0; // core LB method
    virtual void calcMoments(Lattice &lattice) const = 0;
    float getTimeTakenByCollideAndStream() const;
    float getTimeTakenByCalcMoments() const;
    float getTotalTimeTaken() const;

};

#endif
