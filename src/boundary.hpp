#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "pointers.hpp"
#include "agglomerate.hpp"
#include "tuples.hpp"

class Boundary: protected Pointers{
public:
  // 6 entries are for E, W, N, S, U, D bdries
  string type[6];               // can be n - noslip (solid)
                                //        f - freeslip (fluid - not implemented)
                                //        i - inflow   (fluid)
                                //  	  o - outflow  (fluid)
				//        p - periodic
  bool xperiodic, yperiodic, zperiodic;

  string type_th[6];
  double uE[3], uW[3], uN[3];   // bdry velocities
  double uS[3], uU[3], uD[3];

  Boundary(PLBPD *);
  ~Boundary();
  void setup();
  void setup_heatxfer();
  void apply_fluidBC();
  void apply_thermalBC();
  
  // for heat transfer, list of thermal boundary nodes
  // E W N S U D -> 0 1 2 3 4 5
  vector<tuples3<int> > dirichlet[6], neumann[6];

private:
  void calc_ub(int,int,int,int,double *,double *,Agglomerate *,int);

  // local variables to store hydrodynamic force
  double *locforce, *loctorque, *buffer;
};

#endif
