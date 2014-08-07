#ifndef MRT_HPP
#define MRT_HPP

#include "lb_dynamics.hpp"
#include "MersenneTwister/MersenneTwister.h"

class MRT: public LB_Dynamics{
public:
  MRT(PLBPD *);
  ~MRT();
  void setup();
  void init();
  void calc_moments();
  void collideNstream();
  double calc_fluidrhoavg();
  void calc_eqdist(double rholoc, double *uloc, double *nloc);

private:
  double *wsqrt, *winvsqrt;

  // eqlb population distribution
  double *neq;

  // thermal energy
  double kB, temperature;

  // transformation matrix and the vector of the normalization factors
  double **E;

  // normalized functions and their prefactors
  double *x;

  // modes
  double *mode;

  // MRT parameters
  double *gamma, *phi;

  // random numbers
  void gen_random();
  MTRand::uint32 seed[MTRand::N];
  MTRand mt;

  // collideNstream functions
  void cns_no_heatxfer();
  void cns_heat_xfer();
};

#endif
