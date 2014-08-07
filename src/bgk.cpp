#include "bgk.hpp"
#include "domain.hpp"
#include "lattice.hpp"
#include "boundary.hpp"
#include "timer.hpp"
#include "error.hpp"
#include "input.hpp"
#include "model.hpp"
#include "th_dynamics.hpp"

#include <mpi.h>
#include <cmath>
#include <limits>

#include <global.hpp>
using namespace PLBPD_NS; // for machine epsilon

#include "utils.hpp"

/*! \brief description */

/******************************************************************************/

BGK::BGK(PLBPD *plbpd): LB_Dynamics(plbpd){
}

/******************************************************************************/

BGK::~BGK(){
}

/******************************************************************************/

void BGK::setup(){
  nu = input->nu;

  nV = model->fluid->nVelocity;
  c = model->fluid->c;
  w = model->fluid->w;
  cs2 = model->fluid->cs2;

  // set BGK parameter
  omega = 1./(3.*nu+0.5);
  tau = 3.*nu + 0.5;
}

/******************************************************************************/

void BGK::init(){

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  double ****n = lattice->n;
  double ***rho = lattice->rho;
  int ***nodetype = lattice->type;

  double rholoc;
  double *uloc = NULL;
  double *nloc = NULL;

  // initialize dists at domain nodes
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	nloc = n[zl][yl][xl];
	// interior fluid
	if (nodetype[zl][yl][xl]==INTFLD){
	  rholoc = input->rhoFluid;
	  uloc = input->flowVel;
	  calc_eqdist(rholoc, uloc, nloc);
	}
	else if (nodetype[zl][yl][xl]==INFLOW){
	  // inflow nodes (west boundary)
	  rholoc = input->rhoFluid;
	  uloc = boundary->uW;
	  calc_eqdist(rholoc, uloc, nloc);
	}
	else if (nodetype[zl][yl][xl]==OUTFLOW){
	  // outflow nodes (east boundary)
	  rholoc = input->rhoFluid;
	  uloc = boundary->uE;
	  calc_eqdist(rholoc, uloc, nloc);
	}
	else if (nodetype[zl][yl][xl]==NOSLIP){
	  // domain noslip solid
	  // velocities at these nodes have alread been set in
	  // boundary->setup(). do not need to compute nk's
	  rho[zl][yl][xl]=input->rhoSD;
	}
	else if (nodetype[zl][yl][xl]==FREESLIP){
	  string str("'freeslip' boundary condn not implemented");
	  error->all(str);
	}
      } // end for
    } // end for
  } // end for
}

/******************************************************************************/

void BGK::collideNstream(){

  timer->start(COLLIDENSTREAM_FL);
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  double ****n = lattice->n;
  double ****ntmp = lattice->ntmp;
  double ****u = lattice->u;
  double ***rho = lattice->rho;

  double *extF = input->extF;
  double rholoc, ueq[3], usq, neq, cu;

  // non-isothermal case
  if(input->heat_xfer){
    double ***theta = lattice->theta;
    double ***alpha = lattice->alpha;

    double buoF[] = {0.0, 0.0, 0.0};

    double Ra = th_dynamics->Ra;
    double Lx = static_cast<double>(domain->Lx);
    double gbeta_by_alpha = Ra*nu/(Lx*Lx*Lx);
    double theta_film = 0.5;
   
    // at domain and bdry layer nodes
    for(int zl=1; zl<zdim+3; zl++){
      for(int yl=1; yl<ydim+3; yl++){
	for(int xl=1; xl<xdim+3; xl++){
	  rholoc = rho[zl][yl][xl];
	  buoF[2] = rholoc*gbeta_by_alpha*alpha[zl][yl][xl]*
	    (theta[zl][yl][xl]-theta_film);
	  for(int i=0; i<3; i++)
	    ueq[i] = u[zl][yl][xl][i] + (extF[i] + buoF[i])*tau/rholoc;
	  usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
	  for (int k=0; k<nV; k++){
	    // collide and stream
	    cu  = c[k][0]*ueq[0] + c[k][1]*ueq[1] + c[k][2]*ueq[2];
	    neq = w[k]*rholoc*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
	    ntmp[zl+c[k][2]][yl+c[k][1]][xl+c[k][0]][k] =
	      (1.0-omega)*n[zl][yl][xl][k] + omega*neq;
	  }
	}
      }
    }
  }
  // isothermal case
  else{
    // at domain and bdry layer nodes
    for(int zl=1; zl<zdim+3; zl++){
      for(int yl=1; yl<ydim+3; yl++){
	for(int xl=1; xl<xdim+3; xl++){
	  rholoc = rho[zl][yl][xl];
	  for(int i=0; i<3; i++)
	    ueq[i] = u[zl][yl][xl][i] + extF[i]*tau/rholoc;
	  usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
	  for (int k=0; k<nV; k++){
	    // collide and stream
	    cu  = c[k][0]*ueq[0] + c[k][1]*ueq[1] + c[k][2]*ueq[2];
	    neq = w[k]*rholoc*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
	    ntmp[zl+c[k][2]][yl+c[k][1]][xl+c[k][0]][k] =
	      (1.0-omega)*n[zl][yl][xl][k] + omega*neq;
	  }
	}
      }
    }
  }
  
  
  // swap n and ntmp
  lattice->ntmp = n;
  lattice->n = ntmp;
  timer->end(COLLIDENSTREAM_FL);
}
/******************************************************************************/

void BGK::calc_moments(){
  timer->start(CALCMOMENTS_FL);
  // at fluid  (interior, inflow, outflow) nodes

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  double **** const n = lattice->n;
  double **** const u = lattice->u;
  double *** const rho = lattice->rho;
  int *** const nodetype = lattice->type;
  int NODETYPE;

  double uloc[3], *nloc = NULL, rholoc, rhoinv;

  // at domain and boundary fluid nodes
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	NODETYPE = nodetype[zl][yl][xl];

	// fluid nodes
 	if((NODETYPE==INTFLD)or(NODETYPE==INFLOW)or(NODETYPE==OUTFLOW)){
	  nloc = n[zl][yl][xl];
	  rholoc = 0.;
	  uloc[0] = 0.; uloc[1] = 0.; uloc[2] = 0.;
	  for(int k=0; k<nV; k++){
	    rholoc  += nloc[k];
	    uloc[0] += nloc[k]*c[k][0];
	    uloc[1] += nloc[k]*c[k][1];
	    uloc[2] += nloc[k]*c[k][2];
	  }
	  rho[zl][yl][xl] = rholoc;
	  rhoinv = 1./rholoc;
	  for(int i=0; i<3; i++)
	    u[zl][yl][xl][i] = uloc[i]*rhoinv;
	}
      }
    }
  }
  timer->end(CALCMOMENTS_FL);  
}

/******************************************************************************/

double BGK::calc_fluidrhoavg(){
  
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  double ***rho = lattice->rho;
  int ***nodetype = lattice->type;

  double tmp = 0.;
  int ctr = 0;
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	if( (nodetype[zl][yl][xl]==INTFLD) or
	    (nodetype[zl][yl][xl]==INFLOW) or
	    (nodetype[zl][yl][xl]==OUTFLOW) ){
	  ctr++;
	  tmp += rho[zl][yl][xl];
	}
      }
    }
  }

  // sum across all procs
  double rhoTot;
  MPI_Allreduce(&tmp,&rhoTot,1,MPI_DOUBLE,MPI_SUM,world);

  int nfluid;
  MPI_Allreduce(&ctr,&nfluid,1,MPI_INT,MPI_SUM,world);

  return (rhoTot/static_cast<double>(nfluid));
}

/******************************************************************************/

inline void BGK::calc_eqdist(double rholoc, double *uloc, double *nloc){

  double usq, cu, neq;
  usq = uloc[0]*uloc[0] + uloc[1]*uloc[1] + uloc[2]*uloc[2];
  for(int k=0; k<nV; k++){
    cu = c[k][0]*uloc[0] + c[k][1]*uloc[1] + c[k][2]*uloc[2];
    neq = w[k]*rholoc*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    nloc[k] = neq;
  }
}
  
/******************************************************************************/
