// define the model used (e.g. D2Q9, D3Q19 etc.)

#ifndef MODEL_HPP
#define MODEL_HPP

#include "pointers.hpp"

struct Var{
  int nVelocity;
  double cs2;            // square of sound speed
  int (*c)[3];        	 // unit vectors of vel sublattice
  double *w;             // weights for each direction
  //double (*coeff)[4];  // coeffs of (pseudo) eqlb distribution
  int n_unk;             // num of DFs to communicate
  int *unkE, *unkW, *unkN, *unkS, *unkU, *unkD;
  int *reverse;
};

class Model: protected Pointers{
 public:
  Var *fluid;
  Var *thermal;

  // methods
 public:
  Model(PLBPD *);
  ~Model();
  void setup();
};

#endif
