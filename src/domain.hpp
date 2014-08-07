#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "pointers.hpp"

class Domain: protected Pointers{
public:
  // global settings
  int Lx, Ly, Lz;

  // local (my proc) settings
  int lxmin, lxmax, lymin, lymax, lzmin, lzmax;
  int lx, ly, lz;

  Domain(PLBPD *);
  ~Domain();
  void setup();
};

#endif
