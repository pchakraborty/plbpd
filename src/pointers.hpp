#ifndef POINTERS_HPP
#define POINTERS_HPP

#include "mpi.h"
#include "plbpd.hpp"

#include <iostream>
#include <string>
#include <vector>
using std::string;
using std::vector;
using std::cout;
using std::endl;

class Pointers{
public:
  Pointers(PLBPD *ptr):
    plbpd(ptr),
    config(ptr->config),
    input(ptr->input),
    comm(ptr->comm),
    error(ptr->error),
    memory(ptr->memory),
    output(ptr->output),
    restart(ptr->restart),
    timer(ptr->timer),
    model(ptr->model),
    domain(ptr->domain),
    boundary(ptr->boundary),
    lattice(ptr->lattice),
    lb_dynamics(ptr->lb_dynamics),
    th_dynamics(ptr->th_dynamics),
    p_dynamics(ptr->p_dynamics),
    world(ptr->world){}
  virtual ~Pointers(){}

protected:
  PLBPD *plbpd;
  Config *&config;
  Input *&input;
  Comm *&comm;
  Error *&error;
  Memory *&memory;
  Output *&output;
  Restart *&restart;
  Timer *&timer;
  Model *&model;
  Domain *&domain;
  Boundary *&boundary;
  Lattice *&lattice;
  LB_Dynamics *&lb_dynamics;
  TH_Dynamics *&th_dynamics;
  P_Dynamics *& p_dynamics;
  MPI_Comm &world;
};

#endif
  
  
      
