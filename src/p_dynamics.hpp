#ifndef P_DYNAMICS_HPP
#define P_DYNAMICS_HPP

#include "pointers.hpp"
#include "agglomerate.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

class P_Dynamics: protected Pointers{
public:
  unsigned int nAglmrts;      // total number of aglmrts: read from input.plbpd

  P_Dynamics(PLBPD *);
  virtual ~P_Dynamics();
  void init();                // read aglmrt data
  void setup();               // allocate memory
  void integrate();           // time-marching
  Agglomerate *aglmrt;        // structure to hold agglomerate data
                              // all procs maintain this data. since the
                              // num of primaries is not huge, this is OK
  void reset_flag_fluidist(); // node covering/uncovering
  void calc_mass_com_moi();   // calculate the mass, centre of mass and
                              // moment of inertia of the agglomerates
  bool FIRSTCALL2RESET;       // to set node vels according to the aglmrt vels

private:
  void read_aglmrt_data();
  void calc_aglmrtNodeVel(int, int, int, double*, double *, int);
  void checkAglmrtSanity();
  void state2array(Agglomerate *ag, double *y);
  void array2state(Agglomerate *ag, double *y);

  // for ode calc
  int STATE_SIZE;
  const gsl_odeiv_step_type *T;
  gsl_odeiv_step *s;
  double dt; // step size
};

#endif
