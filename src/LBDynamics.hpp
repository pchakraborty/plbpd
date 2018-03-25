#ifndef LBDYNAMICS_HPP
#define LBDYNAMICS_HPP

#include <string>
#include "LBModel.hpp"
#include "Lattice.hpp"

class LBDynamics{

public:

    LBDynamics();
    virtual ~LBDynamics();
    virtual void initialize(Lattice &lattice)=0;
    virtual void collideAndStream(Lattice &lattice)=0; // core LB method
    // virtual void getMoments()=0; // of fluid nodes
    virtual float getAvgFluidDensity()=0;
    
};

#endif
