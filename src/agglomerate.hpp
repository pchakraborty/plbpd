// data structure for agglomerates

#ifndef AGGLOMERATE_HPP
#define AGGLOMERATE_HPP

#include "primary.hpp"
#include <vector>
using std::vector;

class Agglomerate{
public:
  unsigned int ID;           // unique ID
  unsigned int nPrimaries;   // num of primary particles in an agglomerate
  vector<Primary> primary;   // vector containing primary data
  double rho;                // agglomerate density
  double volume;             // agglomerate volume
  double th_diff;            // thermal diffusivity of agglomerates

  /* rigid body variables */
  // constant quantities
  double mass;
  double Ibody[3];           // (pricipal) moments of inertia
  double IbodyInv[3];        // inv(Ibody)

  // state variables
  double com[3];             // centre of mass of agglomerate
  double R[9];               // rotation matrix
  double P[3];               // linear momentum:  m*v
  double L[3];               // angular momentum: I*omega

  // derived quantities (auxiliary variables)
  double Iinv[9];            // inv(Moment of Inertia)
  double u[3];               // linear velocity of centre of mass
  double omega[3];           // angular velocity

  // computed quantities
  double F_hyd[3];           // hydrodynamic force on the agglomerate
  double T_hyd[3];	     // hydrodynamic torque on the agglomerate
  double F_c[3];             // force of covering and uncovering
  double T_c[3];             // torque of covering and uncovering
  double F_tot[3];           // total force
  double T_tot[3];           // total torque
  
  Agglomerate(){
    ID = -1;
    nPrimaries = -1;
    rho = -1.0;
    mass = -1.0;
    for(int i=0; i<3; i++){
      Ibody[i] = 0.0;
      IbodyInv[i] = 0.0;
      com[i] = 0.0;
      P[i] = 0.0;
      L[i] = 0.0;
      u[i] = 0.0;
      omega[i] = 0.0;
      F_hyd[i] = 0.0;
      T_hyd[i] = 0.0;
      F_c[i] = 0.0;
      T_c[i] = 0.0;
      F_tot[i] = 0.0;
      T_tot[i] = 0.0;
    }
    for(int i=0; i<9; i++){
      R[i] = 0.0;
      Iinv[i] = 0.0;
    }
  }
};

#endif
