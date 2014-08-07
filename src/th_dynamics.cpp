#include "th_dynamics.hpp"
#include "domain.hpp"
#include "lattice.hpp"
#include "boundary.hpp"
#include "timer.hpp"
#include "error.hpp"
#include "input.hpp"
#include "model.hpp"
#include "comm.hpp"

#include <mpi.h>
#include <cmath>
#include <limits>
/*! \brief description */

/******************************************************************************/

TH_Dynamics::TH_Dynamics(PLBPD *plbpd): Pointers(plbpd){
}

/******************************************************************************/

TH_Dynamics::~TH_Dynamics(){
}

/******************************************************************************/

void TH_Dynamics::setup(){
  // local copies of model variables
  nV_th = model->thermal->nVelocity;
  
  c = model->thermal->c;
  w = model->thermal->w;
  cs2 = model->thermal->cs2;

  // Rayleigh number
  Ra = input->Ra;
}

/******************************************************************************/

void TH_Dynamics::init(){

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  double **** const n_th = lattice->nth;
  double **** const u = lattice->u;

  double thetaloc;
  double *nloc = NULL, *uloc = NULL;

  double ***alpha = lattice->alpha;

  // initialize dists at domain and boundary nodes
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	thetaloc = 0.;
	// set TH_Dynamics parameter
	alpha[zl][yl][xl] = input->th_diff_fl;

	// u at each node has already been set by
	// lb_dynamics->init()
	uloc = u[zl][yl][xl];

	// compute nth at each node
	nloc = n_th[zl+1][yl+1][xl+1];
	calc_eqdist(thetaloc,uloc,nloc);
      }
    }
  }
}

/******************************************************************************/

void TH_Dynamics::collideNstream(){

  timer->start(COLLIDENSTREAM_TH);
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  double **** const nth = lattice->nth;
  double **** const nth_tmp = lattice->nth_tmp;
  double **** const u = lattice->u;
  double *** const theta = lattice->theta;
  double *** const alpha = lattice->alpha;

  double thetaloc, *uloc = NULL, neq, cu, omega;

  // at domain and bdry layer nodes
  for(int zl=1; zl<zdim+5; zl++){
    for(int yl=1; yl<ydim+5; yl++){
      for(int xl=1; xl<xdim+5; xl++){
	omega = 1.0/(3.0*alpha[zl-1][yl-1][xl-1] + 0.5);
	thetaloc = theta[zl-1][yl-1][xl-1];
	uloc = u[zl-1][yl-1][xl-1];
	for (int k=0; k<nV_th; k++){
	  // collide and stream
	  cu  = c[k][0]*uloc[0] + c[k][1]*uloc[1] + c[k][2]*uloc[2];
	  neq = w[k]*thetaloc*(1.0+3.0*cu);
 	  nth_tmp[zl+c[k][2]][yl+c[k][1]][xl+c[k][0]][k] =
  	    (1.0-omega)*nth[zl][yl][xl][k] + omega*neq;
	}
      }
    }
  }
  
  // swap nth and nth_tmp
  lattice->nth_tmp = nth;
  lattice->nth = nth_tmp;
  timer->end(COLLIDENSTREAM_TH);
}
/******************************************************************************/

void TH_Dynamics::calc_temperature(){
  timer->start(CALCMOMENTS_TH);

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  double *** const theta = lattice->theta;
  double **** const nthloc = lattice->nth;

  for(int zl=1; zl<zdim+3; zl++){
    for(int yl=1; yl<ydim+3; yl++){
      for(int xl=1; xl<xdim+3; xl++){
	double thetaloc = 0.0;
	for(int k=0; k<nV_th; k++)
	  thetaloc += nthloc[zl+1][yl+1][xl+1][k];
	theta[zl][yl][xl] = thetaloc;
      }
    }
  }
  timer->end(CALCMOMENTS_TH);  
}

/******************************************************************************/

double TH_Dynamics::calc_avgtemperature(){
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  double *** const theta = lattice->theta;

  double thetaTot = 0.;
  double nNodes = 0;
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	nNodes += 1.0;
	thetaTot += theta[zl][yl][xl];
      }
    }
  }

  // sum across all procs
  double sendBuf[2]; double recvBuf[2];
  sendBuf[0] = thetaTot; sendBuf[1] = nNodes;
  MPI_Allreduce(&sendBuf,&recvBuf,2,MPI_DOUBLE,MPI_SUM,world);
  return (recvBuf[0]/recvBuf[1]);
}

/******************************************************************************/

void TH_Dynamics::calc_eqdist(double thetaloc, double *uloc, double *nthloc){

  double cu, neq;
  for(int k=0; k<nV_th; k++){
    cu = c[k][0]*uloc[0] + c[k][1]*uloc[1] + c[k][2]*uloc[2];
    neq = w[k]*thetaloc*(1.0+3.0*cu);
    nthloc[k] = neq;
  }
}
  
/******************************************************************************/
