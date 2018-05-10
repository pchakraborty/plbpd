#ifndef LBDYNAMICS_HPP
#define LBDYNAMICS_HPP

#include <string>
#include "LBModel.hpp"
#include "Lattice.hpp"

class LBDynamics{

public:

    LBDynamics();
    virtual ~LBDynamics();
    virtual void collideAndStream(Lattice &lattice)=0; // core LB method
    virtual void calcMoments(Lattice &lattice)=0;
    
};

#endif
