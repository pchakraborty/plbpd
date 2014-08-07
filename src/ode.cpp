#include "agglomerate.hpp"
#include "utils.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

int ode(double t, const double y[], double f[], void *params){
  /* this fn defines the rhs of the ode */

  Agglomerate *ag = (Agglomerate *)params;

  // u
  for(int i=0; i<3; i++)
    *f++ = ag->u[i];

  // dot(R) = (omega*)*(R)
  double Rdot[9], OmegaStar[9];
  utils::star(ag->omega,OmegaStar);
  utils::sqmatmul(OmegaStar,ag->R,Rdot,3);

  for(int i=0; i<9; i++)
    *f++ = Rdot[i];

  for(int i=0; i<3; i++)
    *f++ = ag->F_tot[i];

  for(int i=0; i<3; i++)
    *f++ = ag->T_tot[i];
  
  return GSL_SUCCESS;
}

