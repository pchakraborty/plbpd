#ifndef BGK_HPP
#define BGK_HPP

#include "lb_dynamics.hpp"

class BGK: public LB_Dynamics{
public:
  double nu;

  BGK(PLBPD*);
  ~BGK();
  void setup();
  void init();
  void calc_moments();
  void collideNstream();
  double calc_fluidrhoavg();
  void calc_eqdist(double rholoc, double *uloc, double *nloc);
private:
  double omega;  // parameter for BGK dynamics
  double tau;
};

#endif
