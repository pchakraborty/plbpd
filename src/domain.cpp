#include "domain.hpp"
#include "comm.hpp"
#include "input.hpp"
#include "comm.hpp"

#include "mpi.h"

/******************************************************************************/

Domain::Domain(PLBPD *plbpd): Pointers(plbpd){
  Lx = 0;
  Ly = 0;
  Lz = 0;

  lx = ly = lz = 0;
  lxmin = lxmax = 0;
  lymin = lymax = 0;
  lzmin = lzmax = 0;
}

/******************************************************************************/

Domain::~Domain(){
}

/******************************************************************************/

void Domain::setup(){
  Lx = input->Lx;
  Ly = input->Ly;
  Lz = input->Lz;

  /* divide the global domain into sub-domains */
  lxmin = (comm->myloc[0]  )*Lx/comm->npx;
  lxmax = (comm->myloc[0]+1)*Lx/comm->npx -1;
  lymin = (comm->myloc[1]  )*Ly/comm->npy;
  lymax = (comm->myloc[1]+1)*Ly/comm->npy -1;
  lzmin = (comm->myloc[2]  )*Lz/comm->npz;
  lzmax = (comm->myloc[2]+1)*Lz/comm->npz -1;

  lx = lxmax-lxmin+1;
  ly = lymax-lymin+1;
  lz = lzmax-lzmin+1;
  
  /*
  MPI_Barrier(world);
  if(comm->root){
    std::cout<<
      "proc  E   W   N   S   U   D  xmin xmax ymin ymax zmin zmax lx ly lz\n";
  }
  MPI_Barrier(world);
  printf("%3d %3d %3d %3d %3d %3d %3d %4d %4d %4d %4d %4d %4d %3d %2d %2d\n",
	 comm->me, comm->east, comm->west, comm->north, comm->south, 
	 comm->up, comm->down, lxmin, lxmax, lymin, lymax, lzmin, lzmax,
	 lx, ly, lz);
  MPI_Barrier(world);
  if(comm->root)
    printf("\n");
  */
}

/******************************************************************************/
