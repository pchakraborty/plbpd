#ifndef TH_DYNAMICS_HPP
#define TH_DYNAMICS_HPP

#include "pointers.hpp"

class TH_Dynamics: protected Pointers{
public:
  double Ra;                 // Rayleigh number

  TH_Dynamics(PLBPD *);
  ~TH_Dynamics();

  void setup();
  void init();               // initialize fluid
  void collideNstream();     // core LB method
  void calc_temperature();       
  double calc_avgtemperature(); // compute avg(rho) of fluid
  void calc_eqdist(double thetaloc, double *uloc, double *nloc);

protected:
  // local copies of model variables
  int nV_th;
  int (*c)[3];
  double *w, cs2;
};

#endif
