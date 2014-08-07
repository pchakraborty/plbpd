#include "plbpd.hpp"
#include "mpi.h"
#include <iostream>

int main(int argc, char* argv[]){
  /*! \brief main driver
    @mainpage Parallel Lattice Boltzmann Particle Dynamics (plbpd)
    @author Purnendu Chakraborty, University of Maryland
   */
  MPI_Init(&argc, &argv);

  PLBPD *plbpd;
  try{
    plbpd = new PLBPD(argc,argv,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  catch(std::bad_alloc &xa){
    std::cerr<<"main() failed to allocate memory for plbpd\n";
    return 1;
  }
  
  plbpd->setup();
  plbpd->run();
  delete plbpd;

  MPI_Finalize();

  return 0;
}
