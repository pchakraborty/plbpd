#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "pointers.hpp"

/* node types */
enum{BUFFER=-1, INTFLD=-2, INFLOW=-3, OUTFLOW=-4, NOSLIP=-5, FREESLIP=-6};
// and aglmrt type = i where i is the aglmrt id (i>=0)

/* domain bdry node type: Dirichlet/Neumann (for heat transfer) */
// enum{neumannE, neumannW, neumannN, neumannS, neumannU, neumannD, 
//      dirichletE, dirichletW, dirichletN, dirichletS, dirichletU, dirichletD};

class Lattice: protected Pointers{
public:
  int ***type;      // type of each node
  double ***rho;    // rho at each node
  double ****u;     // u[3] at each node
  double ****n;     // n[nVelocity] at each node [FLUID]
  double ****ntmp;  // ntmp[nVelocity] at each node, for collideNstream()

  // thermal flow
  double ***theta;  // non-dimensional temperature
  double ***alpha;  // thermal diffusivity at each node
  double ****nth;
  double ****nth_tmp;
  double ****RV;    // Random Variables for mrt dynamics

  Lattice(PLBPD *);
  ~Lattice();
  void setup();
};

#endif

