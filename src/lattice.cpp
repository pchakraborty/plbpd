#include "lattice.hpp"
#include "domain.hpp"
#include "boundary.hpp"
#include "memory.hpp"
#include "error.hpp"
#include "model.hpp"
#include "comm.hpp"
#include "input.hpp"

#include "mpi.h"

/******************************************************************************/

Lattice::Lattice(PLBPD *plbpd): Pointers(plbpd){
  type = NULL;
  rho = NULL;
  u = NULL;
  n = NULL;
  ntmp = NULL;
  theta = NULL;
  alpha = NULL;
  nth = NULL;
  nth_tmp = NULL;
  RV = NULL;
}

/******************************************************************************/

Lattice::~Lattice(){
  if(type) memory->delete3Darray(type);
  if(rho) memory->delete3Darray(rho);
  if(u) memory->delete4Darray(u);
  if(n) memory->delete4Darray(n);
  if(ntmp) memory->delete4Darray(ntmp);
  if(theta) memory->delete3Darray(theta);
  if(alpha) memory->delete3Darray(alpha);
  if(nth) memory->delete4Darray(nth);
  if(nth_tmp) memory->delete4Darray(nth_tmp);
  if(RV) memory->delete4Darray(RV);
}

/******************************************************************************/
  
void Lattice::setup(){
  /* create interior fluid and domain bdry (fluid/solid) nodes
     the domain bdry nodes do not need to be set anymore */

  // needs to be called 'after' domain->setup()
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  // each proc builds it own local lattice
  // in addition to the nodes belonging to a given proc, the proc
  // also contains the bdry nodes of the adjacent procs. in addition to that
  // 'n' and 'ntmp' have a buffer layer on each side. that makes their
  // dimensions (zdim+4,ydim+4,xdim+4). type, rho, u can very well
  // have dimensions of (zdim+2,ydim+2,xdim+2), but that makes the
  // coding difficult. for uniformity of indexing, all the arrays have
  // the same dimensions. so,
  // 0:xdim+3 includes both the bdry and buffer layers
  // 1:xdim+2 includes the bdry layer but excludes the buffer layer
  // 2:xdim+1 excludes both the bdry and buffer layers
  // for heat transfer, we need nth and nth_tmp to be (zdim+6,ydim+6,xdim+6)
  // that is because, we need correct value of theta at the bdry layer nodes
  // and not just at the domain nodes.

  type    = memory->create3Darray<int>(zdim+4,ydim+4,xdim+4);
  rho     = memory->create3Darray<double>(zdim+4,ydim+4,xdim+4);
  u       = memory->create4Darray<double>(zdim+4,ydim+4,xdim+4,3);

  int nV_fl = model->fluid->nVelocity;
  n       = memory->create4Darray<double>(zdim+4,ydim+4,xdim+4,nV_fl);
  ntmp    = memory->create4Darray<double>(zdim+4,ydim+4,xdim+4,nV_fl);

  if(input->lb_model=="mrt")
    RV = memory->create4Darray<double>(zdim+4,ydim+4,xdim+4,nV_fl);

  int nV_th = 0;
  if(input->heat_xfer){
    nV_th = model->thermal->nVelocity;
    theta    = memory->create3Darray<double>(zdim+4,ydim+4,xdim+4);
    alpha    = memory->create3Darray<double>(zdim+4,ydim+4,xdim+4);
    nth      = memory->create4Darray<double>(zdim+6,ydim+6,xdim+6,nV_th);
    nth_tmp  = memory->create4Darray<double>(zdim+6,ydim+6,xdim+6,nV_th);
  }

  // init all (domain) nodes to INTFLD
  // and the rest to BUFFER
  for(int zl=0; zl<zdim+4; zl++){
    for(int yl=0; yl<ydim+4; yl++){
      for(int xl=0; xl<xdim+4; xl++){
	rho[zl][yl][xl] = 1.0; // so that there is no div by zero

	// type
	if( (xl>1) and (xl<xdim+2) and
	    (yl>1) and (yl<ydim+2) and
	    (zl>1) and (zl<zdim+2) )
	  type[zl][yl][xl] = INTFLD;
	else
	  type[zl][yl][xl] = BUFFER;

	// u
	for(int i=0; i<3; i++)
	  u[zl][yl][xl][i] = 0.;

	// n
	for(int k=0; k<nV_fl; k++){
	  n[zl][yl][xl][k] = 0.;
	  ntmp[zl][yl][xl][k] = 0.;
	}

	// RV
	if(input->lb_model=="mrt"){
	  for(int i=0; i<nV_fl; i++)
	    RV[zl][yl][xl][i] = 0.0;
	}

	// theta, alpha
	if(input->heat_xfer){
	  theta[zl][yl][xl] = 0.0;
	  alpha[zl][yl][xl] = 0.0;
	}
      }
    }
  }
  
  if(input->heat_xfer){
    for(int zl=0; zl<zdim+6; zl++){
      for(int yl=0; yl<ydim+6; yl++){
	for(int xl=0; xl<xdim+6; xl++){
	  for(int k=0; k<nV_th; k++){
	    nth[zl][yl][xl][k] = 0.;
	    nth_tmp[zl][yl][xl][k] = 0.;
	  }
	}
      }
    }
  }
}

/******************************************************************************/


