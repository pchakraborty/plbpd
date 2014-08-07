#ifndef LB_DYNAMICS_HPP
#define LB_DYNAMICS_HPP

#include "pointers.hpp"

class LB_Dynamics: protected Pointers{
public:
  LB_Dynamics(PLBPD *);
  virtual ~LB_Dynamics();
  virtual void setup()=0;

  virtual void init()=0;               // initialize fluid
  virtual void collideNstream()=0;     // core LB method
  virtual void calc_moments()=0;       // of fluid nodes, return avg(rho)
  virtual double calc_fluidrhoavg()=0; // compute avg(rho) of fluid
  virtual void calc_eqdist(double rholoc, double *uloc, double *nloc)=0;

protected:
  // local copies of model variables
  int nV;
  int (*c)[3];
  double *w, cs2;
};

#endif

  
