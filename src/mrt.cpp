#include "mrt.hpp"
#include "domain.hpp"
#include "lattice.hpp"
#include "boundary.hpp"
#include "timer.hpp"
#include "error.hpp"
#include "input.hpp"
#include "comm.hpp"
#include "model.hpp"
#include "global.hpp"
#include "utils.hpp"
#include "th_dynamics.hpp"
#include "memory.hpp"
#include "plbpd.hpp"

#include <mpi.h>
#include <limits> // for machine epsilon
/******************************************************************************/

MRT::MRT(PLBPD *plbpd): LB_Dynamics(plbpd){
  winvsqrt = NULL;
  wsqrt = NULL;

  neq = NULL;

  E = NULL;
  x = NULL;
  mode = NULL;

  gamma = NULL;
  phi = NULL;
}

/******************************************************************************/

MRT::~MRT(){
  if(wsqrt) delete [] wsqrt;
  if(winvsqrt) delete [] winvsqrt;

  if(neq) delete [] neq;

  if(E) memory->delete2Darray(E);
  if(mode) delete [] mode;
  if(x) delete [] x;

  if(gamma) delete [] gamma;
  if(phi) delete [] phi;
}

/******************************************************************************/

void MRT::setup(){
  kB = input->kB;
  temperature = input->temperature;

  // model variables
  nV = model->fluid->nVelocity;
  c = model->fluid->c;
  w = model->fluid->w;
  cs2 = model->fluid->cs2;
  double onethird = 1.0/3.0;
  if( (cs2>onethird+100.*PLBPD_NS::machineEpsilonD) or
      (cs2<onethird-100.*PLBPD_NS::machineEpsilonD) ){
    string str = "mrt.cpp has been (semi)optimized for cs2=1/3.";
    str += " in the present case, cs2 = " + utils::T2str<double>(cs2);
    error->all(str);
  }

  wsqrt = new double[nV];
  for(int i=0; i<nV; i++)
    wsqrt[i] = sqrt(w[i]);

  winvsqrt = new double[nV];
  for(int i=0; i<nV; i++)
    winvsqrt[i] = 1./wsqrt[i];

  // eqlb popn dist
  neq = new double [nV];

  // build (weighted) transformation matrix E - Adhikari 2005
  E = memory->create2Darray<double>(nV,nV);
  for(int i=0; i<nV; i++){
    double ci_sq = c[i][0]*c[i][0] + c[i][1]*c[i][1] + c[i][2]*c[i][2];

    // k=0:3 correspond to mass and momentum
    E[0][i] = 1.;
    E[1][i] = c[i][0];
    E[2][i] = c[i][1];
    E[3][i] = c[i][2];

    // k=4 correspond to bulk mode. k=5:9 correspond to shear modes
    E[4][i] = ci_sq - 1.;
    E[5][i] = 3.*c[i][0]*c[i][0] - ci_sq;
    E[6][i] = c[i][1]*c[i][1] - c[i][2]*c[i][2];
    E[7][i] = c[i][0]*c[i][1];
    E[8][i] = c[i][1]*c[i][2];
    E[9][i] = c[i][0]*c[i][2];
    
    // k=10:18 correspond to the kinetic (ghost) modes
    E[10][i] = (3.*ci_sq-5.)*c[i][0];
    E[11][i] = (3.*ci_sq-5.)*c[i][1];
    E[12][i] = (3.*ci_sq-5.)*c[i][2];
    E[13][i] = (c[i][1]*c[i][1]-c[i][2]*c[i][2])*c[i][0];
    E[14][i] = (c[i][2]*c[i][2]-c[i][0]*c[i][0])*c[i][1];
    E[15][i] = (c[i][0]*c[i][0]-c[i][1]*c[i][1])*c[i][2];
    E[16][i] = 3.*ci_sq*ci_sq-6.*ci_sq+1.;
    E[17][i] = (2.*ci_sq-3.)*(3.*c[i][0]*c[i][0]-ci_sq);
    E[18][i] = (2.*ci_sq-3.)*(c[i][1]*c[i][1]-c[i][2]*c[i][2]);
  }
  // normalize eigenvectors
  double norm;
  for(int k=0; k<nV; k++){
    // find norm
    norm = 0.;
    for(int i=0; i<nV; i++)
      norm += w[i]*E[k][i]*E[k][i];
    // normalize
    for(int i=0; i<nV; i++)
      E[k][i] *= sqrt(w[i]/norm);
  }

  x = new double[nV];
  mode = new double[nV];

  // // for random number generation
  int inputseed = input->fluct_seed;
  if(inputseed<0){
    string str = "seed for rng must be a non-negative integer";
    error->all(str);
  }
  // each processor generates its unique seed
  // seed with 624*32=19968 bits
  for(int n=0; n<MTRand::N; n++)
    seed[n] = inputseed*n + comm->me;
  mt.seed(seed);

  // set MRT parameters
  gamma = new double[nV];
  phi = new double[nV];
  double lambda_s = -2.0/(6.*input->nu+1);
  double lambda_b = -2.0/(9.*input->nu_b+1);
  double gamma_s = 1. + lambda_s;
  double gamma_b = 1. + lambda_b;
  if((fabs(gamma_s)>=1.) or fabs(gamma_b)>=1.){
    string str = "phi is imaginary! check input values of nu and nu_b";
    error->all(str);
  }
  for(int k=0; k<nV; k++){
    if(k<4){
      gamma[k] = 1.;         // conserved modes
      phi[k] = 0.;
    }
    else if(k==4){
      gamma[k] = gamma_b;    // bulk mode
      phi[k] = sqrt(1.-gamma_b*gamma_b);
    }
    else if((k>4) and (k<10)){
      gamma[k] = gamma_s;    // shear modes
      phi[k] = sqrt(1.-gamma_s*gamma_s);
    }
    else{
      gamma[k] = 0.;         // kinetic (ghost) modes. more accurate bdry
                             // condns can be obtd by tuning these ghost modes
      phi[k] = sqrt(1.-gamma[k]*gamma[k]);
    }
  }
}

/******************************************************************************/

void MRT::init(){
  /*
    initialize n_k at domain nodes
    Eqn 22 of
    Journal of Statistical Physics, Vol. 104, Nos. 5/6, September 2001
      Lattice-Boltzmann Simulations of Particle-Fluid Suspensions
      A. J. C. Ladd and R. Verberg
    OR
    Eqn 4 of
    J. Chem. Phys. 122, 094902 (2005)
      Lattice-Boltzmann dynamics of polymer solutions
      O. B. Usta, A. J. C. Ladd, J. E. Butler
  */

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  double ****n = lattice->n;
  double ***rho = lattice->rho;
  int ***nodetype = lattice->type;

  double rholoc, *nloc = NULL, *uloc = NULL;

  // initialise dists at domain nodes
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

void MRT::calc_moments(){
  timer->start(CALCMOMENTS_FL);
  // at fluid (interior, inflow, outflow) nodes

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  double **** const n = lattice->n;
  double **** const u = lattice->u;
  double *** const rho = lattice->rho;
  int *** const nodetype = lattice->type;
  int NODETYPE;
  double *extF = input->extF;

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
  	    rholoc   += nloc[k];
  	    // Eqn 9 of [2] says to add 0.5*extF[i] to
  	    // rho*uloc[i]; but if i do that the max vel of
  	    // poiseuille flow is off by a factor of 10
  	    uloc[0]  += nloc[k]*c[k][0];
  	    uloc[1]  += nloc[k]*c[k][1];
  	    uloc[2]  += nloc[k]*c[k][2];
  	  }

  	  rho[zl][yl][xl] = rholoc;
  	  rhoinv = 1./rholoc;
  	  for(int i=0; i<3; i++)
  	    u[zl][yl][xl][i] = uloc[i]*rhoinv;
  	}
      } // end for
    } // end for
  } // end for

  timer->end(CALCMOMENTS_FL);

  // generate random numbers for the next step
  // it is being done here, so that the random numbers can
  // be communicated over to the adjacent procs before 
  // the collideNstream step
  gen_random();
}

/******************************************************************************/

void MRT::collideNstream(){
  
  /*
    See papers:
    [1] Journal of Statistical Physics, Vol. 104, Nos. 5/6, September 2001
        Lattice-Boltzmann Simulations of Particle-Fluid Suspensions
	Ladd, R. Verberg
    Or:
    [2] Journal of Chemical Physics 122, 094902 (2005)
        Lattice-Boltzmann simulations of the dynamics of polymer
	solutions in periodic and confined geometries
	Usta, Ladd, Butler
  */

  timer->start(COLLIDENSTREAM_FL);

  // call the appropriate function
  if(input->heat_xfer)
    cns_heat_xfer();
  else
    cns_no_heatxfer();
  
  timer->end(COLLIDENSTREAM_FL);
}

/******************************************************************************/

double MRT::calc_fluidrhoavg(){
  
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

inline void MRT::calc_eqdist(double rholoc, double *uloc, double *nloc){

  double usq, cu, neq1;
  usq = uloc[0]*uloc[0] + uloc[1]*uloc[1] + uloc[2]*uloc[2];
  for(int k=0; k<nV; k++){
    cu = c[k][0]*uloc[0] + c[k][1]*uloc[1] + c[k][2]*uloc[2];
    neq1 = w[k]*rholoc*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    nloc[k] = neq1;
  }
}
  
/******************************************************************************/

void MRT::cns_no_heatxfer(){
  
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  double ****n = lattice->n;
  double ****ntmp = lattice->ntmp;
  double ****u = lattice->u;
  double ***rho = lattice->rho;
  double ****RV = lattice->RV;

  double *extF = input->extF;

  double rholoc, *uloc=NULL, *nloc=NULL, *r=NULL;
  double cu, usq, sqrt_rho, sum, f_i, phi_k;

  double sqrt_mu = sqrt(3.*kB*temperature);

  // at domain and bdry layer nodes
  for(int zl=1; zl<zdim+3; zl++){
    for(int yl=1; yl<ydim+3; yl++){
      for(int xl=1; xl<xdim+3; xl++){
  	// STARTS HERE
  	rholoc = rho[zl][yl][xl];
  	sqrt_rho = sqrt(rholoc);

  	uloc = u[zl][yl][xl];
  	usq = uloc[0]*uloc[0] + uloc[1]*uloc[1] + uloc[2]*uloc[2];
  	for(int i=0; i<nV; i++){
  	  cu = c[i][0]*uloc[0] + c[i][1]*uloc[1] + c[i][2]*uloc[2];
  	  neq[i] = w[i]*rholoc*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
  	}

  	// normalized functions
  	nloc = n[zl][yl][xl];
  	for(int i=0; i<nV; i++)
  	  x[i] = winvsqrt[i]*(nloc[i]-neq[i]);
	
  	// transformation to mode space
  	for(int k=0; k<nV; k++){
  	  sum = 0.;
  	  for(int i=0; i<nV; i++)
  	    sum += E[k][i]*x[i];
  	  mode[k] = sum;
  	}

  	// relaxation
  	r = RV[zl][yl][xl];
  	for(int k=0; k<4; k++)
  	  mode[k] *= gamma[k];
  	for(int k=4; k<nV; k++){
	  phi_k = phi[k]*sqrt_mu*sqrt_rho;
  	  mode[k] = gamma[k]*mode[k] + phi_k*r[k];
	}

  	// back transformation
  	for(int i=0; i<nV; i++){
  	  sum = 0;
  	  for(int k=0; k<nV; k++)
  	    sum += E[k][i]*mode[k];
  	  x[i] = sum;
  	}

  	// new n[i] (add contribution from external force as well)
  	for(int i=0; i<nV; i++){
  	  f_i = 3.*w[i]*(c[i][0]*extF[0] + c[i][1]*extF[1] + c[i][2]*extF[2]);
  	  ntmp[zl+c[i][2]][yl+c[i][1]][xl+c[i][0]][i] = 
  	    neq[i] + wsqrt[i]*x[i] + f_i;
  	}
  	// ENDS HERE
      }
    }
  }

  // swap n and ntmp
  lattice->ntmp = n;
  lattice->n = ntmp;
}

/******************************************************************************/

void MRT::cns_heat_xfer(){

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  double ****n = lattice->n;
  double ****ntmp = lattice->ntmp;
  double ****u = lattice->u;
  double ***rho = lattice->rho;
  double ****RV = lattice->RV;
  double ***theta = lattice->theta;
  double ***alpha = lattice->alpha;

  double *extF = input->extF;

  double rholoc, *uloc=NULL, *nloc=NULL, *r=NULL;
  double cu, usq, sqrt_rho, sum, f_i, phi_k;

  double Ra = th_dynamics->Ra;
  double Lx = static_cast<double>(domain->Lx);
  double prefactor = Ra*input->nu/(Lx*Lx*Lx);
  
  double T_node, sqrt_mu;
  double T_hot = input->T_hot;
  double T_cold = input->T_cold;
  double T_avg = 0.5*(T_hot+T_cold);
  double buoF[] = {0.0, 0.0, 0.0};

  // at domain and bdry layer nodes
  for(int zl=1; zl<zdim+3; zl++){
    for(int yl=1; yl<ydim+3; yl++){
      for(int xl=1; xl<xdim+3; xl++){
  	// STARTS HERE
  	rholoc = rho[zl][yl][xl];
  	sqrt_rho = sqrt(rholoc);

	// node temperature
	T_node = theta[zl][yl][xl]*(T_hot-T_cold) + T_cold;

	// buoyancy force
	buoF[2] = rholoc*prefactor*alpha[zl][yl][xl]*(T_node-T_avg);

	// equilibrium distributions
  	uloc = u[zl][yl][xl];
  	usq = uloc[0]*uloc[0] + uloc[1]*uloc[1] + uloc[2]*uloc[2];
  	for(int i=0; i<nV; i++){
  	  cu = c[i][0]*uloc[0] + c[i][1]*uloc[1] + c[i][2]*uloc[2];
  	  neq[i] = w[i]*rholoc*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
  	}

  	// normalized functions
  	nloc = n[zl][yl][xl];
  	for(int i=0; i<nV; i++)
  	  x[i] = winvsqrt[i]*(nloc[i]-neq[i]);
	
  	// transformation to mode space
  	for(int k=0; k<nV; k++){
  	  sum = 0.;
  	  for(int i=0; i<nV; i++)
  	    sum += E[k][i]*x[i];
  	  mode[k] = sum;
  	}

  	// relaxation
  	r = RV[zl][yl][xl];
  	for(int k=0; k<4; k++)
  	  mode[k] *= gamma[k];
  	for(int k=4; k<nV; k++){
	  sqrt_mu = sqrt(3*kB*T_node);
	  phi_k = phi[k]*sqrt_mu*sqrt_rho;
  	  mode[k] = gamma[k]*mode[k] + phi_k*r[k];
	}

  	// back transformation
  	for(int i=0; i<nV; i++){
  	  sum = 0;
  	  for(int k=0; k<nV; k++)
  	    sum += E[k][i]*mode[k];
  	  x[i] = sum;
  	}

  	// new n[i] (add contribution from external force as well)
  	for(int i=0; i<nV; i++){
  	  f_i = 3.*w[i]*
	    (c[i][0]*extF[0] + c[i][1]*extF[1] + c[i][2]*(extF[2]+buoF[2]));
  	  ntmp[zl+c[i][2]][yl+c[i][1]][xl+c[i][0]][i] = 
	    neq[i] + wsqrt[i]*x[i] + f_i;
  	}
  	// ENDS HERE
      }
    }
  }

  // swap n and ntmp
  lattice->ntmp = n;
  lattice->n = ntmp;
}

/******************************************************************************/

void MRT::gen_random(){
  
  if(input->fluctuations){
    // this is part of collideNstream() although it is being
    // called from calc_moments()
    timer->start(COLLIDENSTREAM_FL);
    
    // // seed the rng at each time step
    // int inputseed = input->fluct_seed;
    // if(inputseed<0){
    //   string str = "seed for rng must be a non-negative integer";
    //   error->all(str);
    // }
    // // each processor generates its unique seed
    // // seed with 624*32=19968 bits
    // for(int n=0; n<MTRand::N; n++)
    //   seed[n] = inputseed*n + comm->me + plbpd->tStep;
    // mt.seed(seed);

    int xdim = domain->lx;
    int ydim = domain->ly;
    int zdim = domain->lz;
    
    double *RV = NULL;
    
    double sqrt3 = sqrt(3);
    double sqrt3x2 = 2*sqrt3;

    // at domain nodes
    for(int zl=2; zl<zdim+2; zl++){
      for(int yl=2; yl<ydim+2; yl++){
	for(int xl=2; xl<xdim+2; xl++){
	  RV = lattice->RV[zl][yl][xl];
	  for(int k=4; k<nV; k++){
	    // RV[k] = mt.randNorm(0.,1.);
	    RV[k] = -sqrt3 + mt.rand(sqrt3x2);
	  }
	}
      }
    } // end for zl

    timer->end(COLLIDENSTREAM_FL);
  }
}
	
/******************************************************************************/
  
