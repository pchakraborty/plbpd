#ifndef LBDYNAMICS_HPP
#define LBDYNAMICS_HPP

#include <string>
#include "LBModel.hpp"
#include "Lattice.hpp"

class LBDynamics{

public:

    LBDynamics();
    virtual ~LBDynamics();
    virtual void setup()=0;
    virtual void collideAndStream()=0; // core LB method
    // virtual void getMoments()=0; // of fluid nodes
    virtual double getAvgFluidDensity()=0;
    virtual void getEqlbDist(double rholoc, const double &uloc, double &nloc)=0;
    
};

#endif
