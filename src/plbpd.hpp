// plbpd.hpp: pointers to fundamental plbpd classes

#ifndef PLBPD_HPP
#define PLBPD_HPP

#include "mpi.h"
#include <string>

class PLBPD{
public:
  class Config *config;
  class Input  *input;             // input script proecessing
  class Comm   *comm;              // inter-processor communication
  class Error  *error;             // error handling
  class Memory *memory;            // memory allocation and memory used
  class Output *output;            // program output
  class Restart *restart;	   // restart read/write
  class Timer  *timer;             // keep track of time
  class Model *model;              // D3Q19, D3Q27 etc.
  class Domain *domain;
  class Boundary *boundary;
  class Lattice *lattice;
  class LB_Dynamics *lb_dynamics;
  class TH_Dynamics *th_dynamics;
  class P_Dynamics *p_dynamics;

  MPI_Comm world;		   // MPI communicator
  MPI_Info info;
  
  PLBPD(int,char **,MPI_Comm,MPI_Info);
  ~PLBPD();
  void destroy();
  void setup();
  void run();

  bool RSTRT;                      // true if restarting a simulation
  std::string infile;

  int timestep;
  int tStep;                       // timestep counter
};

#endif
