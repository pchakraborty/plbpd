#include "p_dynamics.hpp"
#include "lb_dynamics.hpp"
#include "primary.hpp"
#include "agglomerate.hpp"
#include "input.hpp"
#include "lattice.hpp"
#include "error.hpp"
#include "domain.hpp"
#include "timer.hpp"
#include "boundary.hpp"
#include "model.hpp"
#include "utils.hpp"
#include "comm.hpp"

#include "mpi.h"
#include <cstdlib>   // for atof
#include <algorithm> // for min(x,y)
#include <cmath>     // for fabs(x), ceil(x), floor(x)
#include <fstream>
#include <sstream>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include <global.hpp>
using namespace PLBPD_NS; // for machine epsilon

#include "ode.hpp"

using std::min;

/******************************************************************************/

P_Dynamics::P_Dynamics(PLBPD *plbpd): Pointers(plbpd){
  nAglmrts = 0;
  aglmrt = NULL;
  s = NULL;
  T = NULL;
  FIRSTCALL2RESET = true;
  dt = 0.001;
};

/******************************************************************************/

P_Dynamics::~P_Dynamics(){
  if(input->nAg>0){
    if(s) gsl_odeiv_step_free(s);
    if(aglmrt) delete [] aglmrt;
  }
}

/******************************************************************************/

void P_Dynamics::setup(){
  // called after lattice->setup()
  nAglmrts = input->nAg;
  //dt = 1./input->n_pdsteps;
  
  if(nAglmrts>0){
    try{
      aglmrt = new Agglomerate[nAglmrts];
    }
    catch(std::bad_alloc &xa){
      string str("P_Dynamics::setup() failed to allocate memory");
      error->all(str);
    }
    // ode
    STATE_SIZE = 18;
    if(!input->aglmrt_fixed){
      if(input->integrator=="rk4")
	T = gsl_odeiv_step_rk4;   // Classical Runge-Kutta (4th order)
      else if(input->integrator=="rkf45")
	T = gsl_odeiv_step_rkf45; // Runge-Kutta-Fehlberg 45
      else{
	string str("integrator can be one of rk4/rkf45");
	error->all(str);
      }
      s = gsl_odeiv_step_alloc(T,STATE_SIZE);
    }
  }
}

/******************************************************************************/

void P_Dynamics::init(){
  if(nAglmrts>0) read_aglmrt_data();
}

/******************************************************************************/

void P_Dynamics::read_aglmrt_data(){
  // all procs read the data and store all the agglomerate info
  std::ifstream fin;   // fin behaves like cin
  string infile;

  unsigned int numAgFields = 7;
  if(input->heat_xfer)
    numAgFields = 8;
  const unsigned int numPrFields = 4;

  infile = input->aggDataFile;
  fin.open(infile.c_str());
  if (!fin){
    string str("agglomerate data file " + infile + " could not be opened");
    error->all(str);
  }

  /* read agglomerate data and store them */
  unsigned int prmryID = 0;
  for(unsigned int iAg=0; iAg<nAglmrts; iAg++){
    // read the first 2 lines for aglmrt data and num primaries
    // in the input file, data needs to be in order, i.e. sorted.
    // agglomerate IDs should start with 0 and increase by 1
    string line;
    getline(fin,line);          // velocities (lin, ang), density
    if(line.empty()){
      fin.close();
      string str("aglmrt data (velocity etc.) not found for aglmrt " +
		 utils::T2str<int>(iAg));
      error->all(str);
    }
    else{
      vector<string> words = utils::split(line);
      if(words.size()<numAgFields){
	fin.close();
	string str(utils::T2str<int>(words.size())+"values read (instead of " +
		   utils::T2str<int>(numAgFields) +") on the first line of " +
		   "agglomerate " + utils::T2str<int>(iAg));
	error->all(str);
      }
      // store
      aglmrt[iAg].ID       = iAg;
      aglmrt[iAg].u[0]     = utils::str2T<double>(words[0]);
      aglmrt[iAg].u[1]     = utils::str2T<double>(words[1]);
      aglmrt[iAg].u[2]     = utils::str2T<double>(words[2]);
      aglmrt[iAg].omega[0] = utils::str2T<double>(words[3]);
      aglmrt[iAg].omega[1] = utils::str2T<double>(words[4]);
      aglmrt[iAg].omega[2] = utils::str2T<double>(words[5]);
      aglmrt[iAg].rho      = utils::str2T<double>(words[6]);
      if(input->heat_xfer)
	aglmrt[iAg].th_diff  = utils::str2T<double>(words[7]);
    }

    getline(fin,line);          // 2nd line: num of primaries
    if(line.empty()){
      fin.close();
      string str("num primaries not found for agglomerate " 
		 + utils::T2str<int>(iAg));
      error->all(str);
    }
    else aglmrt[iAg].nPrimaries = atoi(line.c_str());

    if(aglmrt[iAg].nPrimaries<1){
      string str = "aglmrt " + utils::T2str<int>(iAg) + " has "+
	utils::T2str<int>(aglmrt[iAg].nPrimaries) + " primaries";
      error->all(str);
    }
    if(aglmrt[iAg].nPrimaries>PLBPD_NS::MAXPRMRYS){
      string str = "aglmrt " + utils::T2str<int>(iAg) + " has " +
	utils::T2str<int>(aglmrt[iAg].nPrimaries) + "primaries!! You need to" +
	"increase the value of the variable 'MAXPRMRYS' in 'global.hpp' " +
	"from " + utils::T2str<int>(PLBPD_NS::MAXPRMRYS) + " to a value " +
	"that is >= " + utils::T2str<int>(aglmrt[iAg].nPrimaries) +
	". This is important for writing restart files";
      error->all(str);
    }
    
    // next read and store the primary data
    int ctr = 0;
    for(unsigned int iPr=0; iPr<aglmrt[iAg].nPrimaries; iPr++){
      getline(fin,line);
      if(line.empty()){
	fin.close();
	string str="num primaries not found for aglmrt "+utils::T2str<int>(iAg);
	error->all(str);
      }
      ctr++;
      vector<string> words = utils::split(line);
      if(words.size()<numPrFields){
	fin.close();
	string str = utils::T2str<int>(words.size()) + " values read " +
	  "(instead of " + utils::T2str<int>(words.size()) + ") for primary " +
	  utils::T2str<int>(ctr) + " of aglmrt " + utils::T2str<int>(iAg);
	error->all(str);
      }

      // store
      Primary tmp;
      tmp.ID   = prmryID;
      tmp.r    = atof(words[0].c_str());
      tmp.c[0] = atof(words[1].c_str());
      tmp.c[1] = atof(words[2].c_str());
      tmp.c[2] = atof(words[3].c_str());
      
      aglmrt[iAg].primary.push_back(tmp);

      prmryID++;
    }

    // check to see if all primaries (for the given aglmrt) have been read
    if(aglmrt[iAg].primary.size()<aglmrt[iAg].nPrimaries){
      int size = aglmrt[iAg].primary.size();
      fin.close();
      string str = infile + " has " + utils::T2str<int>(size) +
	" primaries (instead of " + utils::T2str<int>(aglmrt[iAg].nPrimaries) +
	") for aglmrt " + utils::T2str<int>(iAg);
      error->all(str);
    }
  }

  fin.close();
}

/******************************************************************************/

void P_Dynamics::reset_flag_fluidist(){
  // all procs do this
  
  timer->start(RESET_FLAG_MOMENT);

  /* 
     first check aglmrt sanity;
     consecutive primaries in the agglomerate must either overlap or touch
     i.e. ||c1c2|| <= r1 + r2 AND ||c1c2|| > r1, ||c1c2|| > r2
  */
  //checkAglmrtSanity();

  double xg, yg, zg;             // global coords of a node
  double dx, dy, dz;             // distance of node from primary centre
  double rP, cP[3];              // local variables for primary data

  int newType;
  int oldType;
  int AGLMRTYPE;

  int nV = model->fluid->nVelocity;
  int (*c)[3] = model->fluid->c;

  int nx, ny, nz;                // nbring node coordinates

  int ***nodetype = lattice->type;
  double ***rho   = lattice->rho;
  double ****u    = lattice->u;
  double ****n    = lattice->n;
  double ***alpha = lattice->alpha;

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  double Lx = static_cast<double>(domain->Lx);
  double Ly = static_cast<double>(domain->Ly);
  double Lz = static_cast<double>(domain->Lz);

  bool xperiodic = boundary->xperiodic;
  bool yperiodic = boundary->yperiodic;
  bool zperiodic = boundary->zperiodic;

  // set forces and torques due to covering/uncovering to zero
  for(int iAg=0; iAg<nAglmrts; iAg++){
    for(int i=0; i<3; i++){
      aglmrt[iAg].F_c[i] = 0.0;
      aglmrt[iAg].T_c[i] = 0.0;
    }
  }

  // sweep the nodes to see which nodes remain aglmrt nodes,
  // switch from intfld to aglmrt node or vice-versa
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	oldType = nodetype[zl][yl][xl];
	newType = INTFLD;
	// global coords of current node, 0<=xg<=Lx-1 etc.
	// xl=2, yl=2, zl=2 in proc 0 has coords (0,0,0)
	xg = static_cast<double>((xl - 2) + domain->lxmin);
	yg = static_cast<double>((yl - 2) + domain->lymin);
	zg = static_cast<double>((zl - 2) + domain->lzmin);
	
	/*************************************************/
	// check if (xl,yl,zl) is inside any aglmrt
	for(unsigned int iAg=0; iAg<nAglmrts; iAg++){
	  AGLMRTYPE = aglmrt[iAg].ID;
	  for(unsigned int iPr=0; iPr<aglmrt[iAg].primary.size(); iPr++){
	    rP    = aglmrt[iAg].primary[iPr].r;    // radius
	    cP[0] = aglmrt[iAg].primary[iPr].c[0]; // (global) coords of centre
	    cP[1] = aglmrt[iAg].primary[iPr].c[1];
	    cP[2] = aglmrt[iAg].primary[iPr].c[2];
	    // correction for periodic bdry
	    if(xperiodic) cP[0] = utils::mod(cP[0],Lx);
	    if(yperiodic) cP[1] = utils::mod(cP[1],Ly);
	    if(zperiodic) cP[2] = utils::mod(cP[2],Lz);

	    // distance of the node from primary centre
	    dx = xg-cP[0]; dy = yg-cP[1]; dz = zg-cP[2];
	    // correction for periodic bdry
	    if(xperiodic) dx -= utils::nint(dx/Lx)*Lx;
	    if(yperiodic) dy -= utils::nint(dy/Ly)*Ly;
	    if(zperiodic) dz -= utils::nint(dz/Lz)*Lz;

	    if((dx*dx + dy*dy + dz*dz) < rP*rP){
	      // node covering
	      newType = AGLMRTYPE;
	      nodetype[zl][yl][xl] = newType;

	      // old values of rho and u
	      double oldrho = rho[zl][yl][xl];
	      double oldu[3];
	      for(int i=0; i<3; i++)
		oldu[i] = u[zl][yl][xl][i];

	      // new values of rho and u
	      rho[zl][yl][xl] = aglmrt[iAg].rho;   // set rho
	      double rprime[3], Uprime[3];
	      calc_aglmrtNodeVel(xg,yg,zg,rprime,Uprime,iAg);
	      for(int i=0; i<3; i++)
		u[zl][yl][xl][i] = Uprime[i];      // set u

	      // // force and torque of covering
	      // if(!FIRSTCALL2RESET){
	      // 	// cout<<"node covered\n";
	      // 	// force = old momentum - new momentum
	      // 	double force[3];
	      // 	for(int i=0; i<3; i++){
	      // 	  force[i] = oldrho*oldu[i] - aglmrt[iAg].rho*Uprime[i];
	      // 	  aglmrt[iAg].F_c[i] += force[i];
	      // 	}
	      // 	// torque
	      // 	double torque[3];
	      // 	for(int i=0; i<3; i++){
	      // 	  utils::cross(rprime, force, torque);
	      // 	  aglmrt[iAg].T_c[i] += torque[i];
	      // 	}
	      // }

	      // if this is the first call to reset__flag_fluidist()
	      // set the node velocities to the aglmrt velocities
	      if(FIRSTCALL2RESET)
		for(int i=0; i<3; i++)
		  u[zl][yl][xl][i] = aglmrt[iAg].u[i];

	      // error check
	      if((oldType==INFLOW) or (oldType==OUTFLOW) or
		 (oldType==NOSLIP) or (oldType==FREESLIP)){
		string str = "primary " + utils::T2str<int>(iPr) +
		  " of aglmrt " + utils::T2str<int>(iAg) + " is at " +
		  "INFLOW/OUTFLOW/NOSLIP/FREESLIP";
		error->all(str);
	      }
	      if((oldType>=0) and (oldType!=newType)){
		string str = "aglmrts " + utils::T2str<int>(oldType) +
		  " and " + utils::T2str<int>(newType) +
		  " are too close to each other";
	      }
	      // set thermal diffusivity of agglomerate
	      if(input->heat_xfer)
		alpha[zl][yl][xl] = aglmrt[iAg].th_diff;
	    }
	  } // end iPr
	} // end iAg
	/*************************************************/

	// node uncovering
	if((oldType>=0) and (newType==INTFLD)){
	  nodetype[zl][yl][xl] = INTFLD; 
	  double oldrho = rho[zl][yl][xl];
 
	  // generate nk:
	  // rho is an avrg of nbring INTFLD node densities
	  double tmprho = 0.;
	  double ctr = 0.;
	  for(int k=1; k<nV; k++){ // cannot include the current node, of course
	    nx = xl + c[k][0];
	    ny = yl + c[k][1];
	    nz = zl + c[k][2];
	    if(nodetype[nz][ny][nx]==INTFLD){
	      tmprho += rho[nz][ny][nx];
	      ctr += 1.0;
	    }
	  }
	  if(fabs(ctr)<100.*PLBPD_NS::machineEpsilonD){
	    string str("uncovered node does not have INTFLD nbrs");
	    error->all(str);
	  }
	  tmprho = tmprho/ctr;
	  // since we update the aglmrt node velocities at each time
	  // step, we already know the u at this node
	  lb_dynamics->calc_eqdist(tmprho, u[zl][yl][xl], n[zl][yl][xl]);

	  // // force and torque of uncovering: (oldrho - newrho)*u
	  // // force = old momentum - new momentum
	  // double force[3];
	  // // force
	  // for(int i=0; i<3; i++){
	  //   force[i] = (oldrho - tmprho)*u[zl][yl][xl][i];
	  //   aglmrt[oldType].F_c[i] += force[i];
	  // }
	  // // torque
	  // double torque[3], rprime[3], R[3];
	  // for(int i=0; i<3; i++)
	  //   R[i] = aglmrt[oldType].com[i];
	  // if(xperiodic) R[0] = utils::mod(R[0],Lx);  
	  // if(yperiodic) R[1] = utils::mod(R[1],Ly);
	  // if(zperiodic) R[2] = utils::mod(R[2],Lz);
	  // double tmp[3];
	  // for(int i=0; i<3; i++)
	  //   tmp[i] = rprime[i] - R[i];
	  // if(xperiodic) tmp[0] -= utils::nint(tmp[0]/Lx)*Lx;
	  // if(yperiodic) tmp[1] -= utils::nint(tmp[1]/Ly)*Ly;
	  // if(zperiodic) tmp[2] -= utils::nint(tmp[2]/Lz)*Lz;
	  // for(int i=0; i<3; i++)
	  //   rprime[i] = aglmrt[oldType].com[i] + tmp[i];
	  // for(int i=0; i<3; i++){
	  //   utils::cross(rprime, force, torque);
	  //   aglmrt[oldType].T_c[i] += torque[i];
	  // }

	  // set thermal diffusivity back to that of the fluid
	  if(input->heat_xfer)
	    alpha[zl][yl][xl] = input->th_diff_fl;
	}
      }
    }
  }
  FIRSTCALL2RESET = false;

  timer->end(RESET_FLAG_MOMENT);
}

/******************************************************************************/

void P_Dynamics::integrate(){
  timer->start(INTEGRATE);

  // SOLVE ODE: dy_i/dt = f_i
  double y[STATE_SIZE], y_err[STATE_SIZE];
  double dydt_in[STATE_SIZE], dydt_out[STATE_SIZE];

  for(unsigned int iAg=0; iAg<nAglmrts; iAg++){
    // external force (acting at the com of the aglmrt - so torque
    // about an axis through the com is 0)
    double extFaglmrt[3] = {0.0, 0.0, 0.0};
    for(int i=0; i<3; i++)
      extFaglmrt[i] = aglmrt[iAg].volume*input->extFaglmrt[i];

    // total force and torque on the aglmrt    
    for(int i=0; i<3; i++){
      aglmrt[iAg].F_tot[i] = aglmrt[iAg].F_hyd[i] + extFaglmrt[i];
      aglmrt[iAg].T_tot[i] = aglmrt[iAg].T_hyd[i];
    }

    state2array(&aglmrt[iAg], &y[0]); // the state variables are in array 'y'

    gsl_odeiv_system dydt = {ode, NULL, STATE_SIZE, &aglmrt[iAg]};
    
    double t = 0.0;
    double t_end = 1.0;

    // intialize dydt_in from system parameters
    GSL_ODEIV_FN_EVAL(&dydt,t,y,dydt_in);

    // finally, integrate
    while(t<t_end){
      int status = gsl_odeiv_step_apply(s,t,dt,y,y_err,dydt_in,dydt_out,&dydt);
      if(status!=GSL_SUCCESS){
	string str("ode integration failed");
	error->all(str);
      }
      for(int i=0; i<STATE_SIZE; i++)
	dydt_in[i] = dydt_out[i];
      
      t += dt;
    }
    array2state(&aglmrt[iAg], &y[0]);
  }
  
  timer->end(INTEGRATE);
}

/******************************************************************************/

void P_Dynamics::checkAglmrtSanity(){
  /* consecutive primaries in the agglomerate must either overlap or touch
     i.e. ||c1c2|| <= r1 + r2 AND ||c1c2|| > r1, ||c1c2|| > r2 */

  // TODO: MAKE SURE 2 AGLMRTS ARE NOT CLOSE

  double r1, r2, c1[3], c2[3], c1c2_sq, r1_sq, r2_sq, r1pr2_sq;

  for(unsigned int iAg=0; iAg<nAglmrts; iAg++){
    for(unsigned int iPr=0; iPr<aglmrt[iAg].primary.size()-1; iPr++){
      // can handle aglmrts with 1 primary
      r1    = aglmrt[iAg].primary[iPr].r;    // radius
      c1[0] = aglmrt[iAg].primary[iPr].c[0]; // (global) coords of centre
      c1[1] = aglmrt[iAg].primary[iPr].c[1];
      c1[2] = aglmrt[iAg].primary[iPr].c[2];

      r2    = aglmrt[iAg].primary[iPr+1].r; 
      c2[0] = aglmrt[iAg].primary[iPr+1].c[0];
      c2[1] = aglmrt[iAg].primary[iPr+1].c[1];
      c2[2] = aglmrt[iAg].primary[iPr+1].c[2];
      
      r1_sq = r1*r1;
      r2_sq = r2*r2;
      r1pr2_sq = (r1+r2)*(r1+r2);
      c1c2_sq = 
	(c1[0]-c2[0])*(c1[0]-c2[0]) +
	(c1[1]-c2[1])*(c1[1]-c2[1]) +
	(c1[2]-c2[2])*(c1[2]-c2[2]);

      if((c1c2_sq>r1pr2_sq) or (c1c2_sq<=r1) or (c1c2_sq<=r2)){
	string str = "prmries " + utils::T2str<int>(iPr) + "and " +
	  utils::T2str<int>(iPr+1) + "of aglmrt " + utils::T2str<int>(iAg) +
	  " are either too far from each other or " +
	  "one is contained in the other";
	error->all(str);
      }
    }
  }
}

/******************************************************************************/

void P_Dynamics::calc_aglmrtNodeVel(int xg, int yg, int zg,
				    double *rp, double *up, int iAg){

  double Lx = static_cast<double>(domain->Lx);
  double Ly = static_cast<double>(domain->Ly);
  double Lz = static_cast<double>(domain->Lz);

  // xg, yg, zg are box coordinates
  rp[0] = static_cast<double>(xg);
  rp[1] = static_cast<double>(yg);
  rp[2] = static_cast<double>(zg);

  bool xperiodic = boundary->xperiodic;
  bool yperiodic = boundary->yperiodic;
  bool zperiodic = boundary->zperiodic;

  /* aglmrt node velocity
     up = U + Omega x (rp-X) .. Eqn 2.8, Ladd 1994b
     where:
       U:     lin vel of the com of the aglmrt iAg
       Omega: ang vel about the com of iAg
       X:     com of iAg */
  double X[3], *U=NULL, *Omega=NULL;
  X[0] = aglmrt[iAg].com[0];
  X[1] = aglmrt[iAg].com[1];
  X[2] = aglmrt[iAg].com[2];
  // correction for periodic bdry
  if(xperiodic) X[0] = utils::mod(X[0],Lx);  
  if(yperiodic) X[1] = utils::mod(X[1],Ly);
  if(zperiodic) X[2] = utils::mod(X[2],Lz);

  U = aglmrt[iAg].u;
  Omega = aglmrt[iAg].omega;
  
  // tmp = rp - X
  double tmp[3];
  for(int i=0; i<3; i++)
    tmp[i] = rp[i] - X[i];
  // correction for periodic bdry
  if(xperiodic) tmp[0] -= utils::nint(tmp[0]/Lx)*Lx;
  if(yperiodic) tmp[1] -= utils::nint(tmp[1]/Ly)*Ly;
  if(zperiodic) tmp[2] -= utils::nint(tmp[2]/Lz)*Lz;

  // compute up
  up[0] = U[0] + Omega[1]*tmp[2] - Omega[2]*tmp[1];
  up[1] = U[1] - Omega[0]*tmp[2] + Omega[2]*tmp[0];
  up[2] = U[2] + Omega[0]*tmp[1] - Omega[1]*tmp[0];
}

/******************************************************************************/

void P_Dynamics::calc_mass_com_moi(){
  // this fn is called only once, right after the first call
  // to reset_flag_dist(). it computes the mass and moi (I) of the aglmrts
  // at time 0. mass of a given aglmrt is constant for the rest of
  // the simulation. Also, I(t) = R(t)*Ibody*R(t)^T with Ibody being
  // constant for the rest of the simulation. this fn sets Ibody and the
  // initial value of R, R(0), which is of course the matrix of eigen
  // vectors of I.

  // NOTE: INITIALLY PRIMARIES CANNOT CROSS BOUNDARY

  int ***nodetype = lattice->type;

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  double xg, yg, zg;          // global coords of node

  double *M = new double[nAglmrts];         // mass
  double *Vol = new double[nAglmrts];       // volume
  double *X = new double[3*nAglmrts];       // com
  double *I = new double[9*nAglmrts];       // moi
  double *sbuffer = new double[9*nAglmrts]; // buffer for adding across procs

  // initialize
  for(unsigned int iAg=0; iAg<nAglmrts; iAg++){
    M[iAg] = 0.0;
    Vol[iAg] = 0.0;
    for(int j=0; j<3; j++)
      X[3*iAg+j] = 0.0;
    for(int j=0; j<9; j++){
      I[9*iAg+j] = 0.0;
      sbuffer[9*iAg+j] = 0.0;
    }
  }

  // first compute com, Vol, M
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	int iAg = nodetype[zl][yl][xl];

	if(iAg>=0){ // aglmrt node
	  // global coords of current node, 0<=xg<=Lx-1 etc.
	  // xl=2, yl=2, zl=2 in proc 0 has coords (0,0,0)
	  xg = static_cast<double>((xl-2)+domain->lxmin);
	  yg = static_cast<double>((yl-2)+domain->lymin);
	  zg = static_cast<double>((zl-2)+domain->lzmin);
	  
	  double rho = aglmrt[iAg].rho;
	  M[iAg] += rho;
	  Vol[iAg] += 1.;

	  X[3*iAg+0] += rho*xg;   // volume = 1
	  X[3*iAg+1] += rho*yg;
	  X[3*iAg+2] += rho*zg;
	}
      }
    }
  }
  
  // sum the quantities over all procs
  // we need the info in all the procs, hence Allreduce
  int count = nAglmrts;   // mass
  for(int i=0; i<count; i++)
    sbuffer[i] = M[i];
  MPI_Allreduce(sbuffer,M,count,MPI_DOUBLE,MPI_SUM,world);
  count = nAglmrts;       // volume
  for(int i=0; i<count; i++)
    sbuffer[i] = Vol[i];
  MPI_Allreduce(sbuffer,Vol,count,MPI_DOUBLE,MPI_SUM,world);
  count = 3*nAglmrts;     // centre of mass
  for(int i=0; i<count; i++)
    sbuffer[i] = X[i];
  MPI_Allreduce(sbuffer,X,count,MPI_DOUBLE,MPI_SUM,world);
  for(int iAg=0; iAg<nAglmrts; iAg++)
    for(int j=0; j<3; j++)
      X[3*iAg+j] = X[3*iAg+j]/M[iAg];
  
  // next compute moi
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	int iAg = nodetype[zl][yl][xl];

	if(iAg>=0){ // aglmrt node
	  // global coords of current node, 0<=xg<=Lx-1 etc.
	  // xl=2, yl=2, zl=2 in proc 0 has coords (0,0,0)

	  // CAREFUL: initially aglmrts should not cross periodic bdrys
	  xg = static_cast<double>((xl-2)+domain->lxmin) - X[3*iAg+0];
	  yg = static_cast<double>((yl-2)+domain->lymin) - X[3*iAg+1];
	  zg = static_cast<double>((zl-2)+domain->lzmin) - X[3*iAg+2];

	  // cout<<"("<<xl<<","<<yl<<","<<zl<<"): ";
	  // cout<<"("<<xg<<","<<yg<<","<<zg<<")"<<endl;

	  // correction for periodicity
	  
	  double rho = aglmrt[iAg].rho;

	  double Ixx =  rho*(yg*yg + zg*zg);
	  double Iyy =  rho*(xg*xg + zg*zg);
	  double Izz =  rho*(xg*xg + yg*yg);
	  double Ixy = -rho*xg*yg;
	  double Ixz = -rho*xg*zg;
	  double Iyz = -rho*yg*zg;

	  I[9*iAg+0] += Ixx;
	  I[9*iAg+1] += Ixy;
	  I[9*iAg+2] += Ixz;
	  I[9*iAg+3] += Ixy;
	  I[9*iAg+4] += Iyy;
	  I[9*iAg+5] += Iyz;
	  I[9*iAg+6] += Ixz;
	  I[9*iAg+7] += Iyz;
	  I[9*iAg+8] += Izz;
	}
      }
    }
  }
  // add moi over procs
  count = 9*nAglmrts;     // moment of inertia
  for(int i=0; i<count; i++)
    sbuffer[i] = I[i];
  MPI_Allreduce(sbuffer,I,count,MPI_DOUBLE,MPI_SUM,world);

  // allocate space for eigval and eigvec computaions
  gsl_vector *Lambda = gsl_vector_alloc(3);
  gsl_matrix *V = gsl_matrix_alloc(3,3);
  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(3); // workspace

  for(unsigned int iAg=0; iAg<nAglmrts; iAg++){
    aglmrt[iAg].mass = M[iAg];
    aglmrt[iAg].volume = Vol[iAg];
    for(int i=0; i<3; i++)
      aglmrt[iAg].com[i] = X[3*iAg+i];

    /* Diagonalize the symmetric matrix I, i.e. compute R and Ibody
       such that I = R*Ibody*R^T, where R is the rotation matrix */
    double data[9];  // store I of iAg in data
    for(int i=0; i<9; i++)
      data[i] = I[9*iAg+i];

    // reshape 'data' into a 3x3 gsl_matrix
    gsl_matrix_view Mat = gsl_matrix_view_array(data,3,3);
    
    // compute eigenvalues and eigenvectors
    gsl_eigen_symmv(&Mat.matrix,Lambda,V,w);

    // sort eig vals and eig vecs (ascending order in magnitude)
    gsl_eigen_symmv_sort(Lambda,V,GSL_EIGEN_SORT_ABS_ASC);

    // store Lambda in agglomerate->Ibody
    for(int i=0; i<3; i++){
      aglmrt[iAg].Ibody[i] = gsl_vector_get(Lambda,i);
      aglmrt[iAg].IbodyInv[i] = 1./aglmrt[iAg].Ibody[i];
    }
    
    // store V in agglomerate->R
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
	aglmrt[iAg].R[i*3+j] = gsl_matrix_get(V,i,j);
    
    // RT = R^T
    double RT[9]; utils::transpose(aglmrt[iAg].R,RT,3);

    /*
    // check
    // R*RT = RT*R = I
    double RRT[9], RTR[9];
    utils::sqmatmul(aglmrt[iAg].R,RT,RRT,3);
    utils::sqmatmul(RT,aglmrt[iAg].R,RTR,3);
    if(comm->me==0){
      cout<<"\nR*R^T:\n";
      utils::printSqMatrix(RRT,3);
      cout<<"\nR^T*R:\n";
      utils::printSqMatrix(RTR,3);
    }
    // I = R*Ibody*R^T
    double tmp1[9], tmp2[9];
    // expand IbodyInv into the diag matrix
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
    	if(i==j) tmp1[3*i+j] = aglmrt[iAg].Ibody[i];
    	else tmp1[3*i+j] = 0.0;
      }
    }
    utils::sqmatmul(aglmrt[iAg].R,tmp1,tmp2,3);
    utils::sqmatmul(tmp2,RT,tmp1,3);
    if(comm->me==0){
      cout<<"\nR*Ibody*R^T:\n";
      utils::printSqMatrix(tmp1,3);
    }
    */

    // P (linear momentum) = m*u
    for(int i=0; i<3; i++)
      aglmrt[iAg].P[i] = M[iAg]*aglmrt[iAg].u[i];

    // L (angular momentum) = I*omega (I is already stored in data)
    utils::matvecprd(data,aglmrt[iAg].omega,aglmrt[iAg].L,3);

    // Iinv = R*IbodyInv*R^T
    // expand IbodyInv into the diag matrix
    double diag[9], tmpMatrix[9];
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	if(i==j) diag[3*i+j] = aglmrt[iAg].IbodyInv[i];
	else diag[3*i+j] = 0.0;
      }
    }
    utils::sqmatmul(aglmrt[iAg].R,diag,tmpMatrix,3);  // tmpMatrix = R*IbodyInv
    utils::sqmatmul(tmpMatrix,RT,aglmrt[iAg].Iinv,3); // Iinv = R*IbodyInv*R^T

    // finally compute primary centres in body space
    double tmpVector[3];
    for(unsigned int iPr=0; iPr<aglmrt[iAg].primary.size(); iPr++){
      for(int i=0; i<3; i++)
	tmpVector[i] = aglmrt[iAg].primary[iPr].c[i] - aglmrt[iAg].com[i];
      utils::matvecprd(RT,tmpVector,aglmrt[iAg].primary[iPr].cb,3);
    }
  }

  /* free up space allocations */
  gsl_eigen_symmv_free(w); // workspace
  gsl_vector_free(Lambda);
  gsl_matrix_free(V);
  delete [] M;
  delete [] Vol;
  delete [] X;
  delete [] I;
  delete [] sbuffer;
}

/******************************************************************************/
 
void P_Dynamics::state2array(Agglomerate *ag, double *y){
  // coords of com
  for(int i=0; i<3; i++)
    *y++ = ag->com[i];

  // rotation matrix R
  for(int i=0; i<9; i++)
    *y++ = ag->R[i];
  
  // linear momentum
  for(int i=0; i<3; i++)
    *y++ = ag->P[i];
  
  // angular momentum
  for(int i=0; i<3; i++)
    *y++ = ag->L[i];
}

/******************************************************************************/

void P_Dynamics::array2state(Agglomerate *ag, double *y){
  // com
  for(int i=0; i<3; i++)
    ag->com[i] = *y++; // new com

  // R
  for(int i=0; i<9; i++)
    ag->R[i] = *y++;

  // P
  for(int i=0; i<3; i++)
    ag->P[i] = *y++;

  // L
  for(int i=0; i<3; i++)
    ag->L[i] = *y++;

  // compute auxiliary variables
  // v = P/mass
  double mass = ag->mass;
  for(int i=0; i<3; i++)
    ag->u[i] = ag->P[i]/mass;

  // Iinv = R*IbodyInv*R^T
  double RT[9];
  utils::transpose(ag->R,RT,3);
  double IbodyInv[9];
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      if(i==j)
	IbodyInv[3*i+j] = ag->IbodyInv[i];
      else
	IbodyInv[3*i+j] = 0.0;
    }
  }
  double prd[9];
  utils::sqmatmul(IbodyInv,RT,prd,3);
  utils::sqmatmul(ag->R,prd,ag->Iinv,3);

  // omega = Iinv*L
  utils::matvecprd(ag->Iinv,ag->L,ag->omega,3);

  // finally update the primary centres
  double tmp[3];
  for(unsigned int iPr=0; iPr<ag->primary.size(); iPr++){
    utils::matvecprd(ag->R,ag->primary[iPr].cb,tmp,3);
    for(int i=0; i<3; i++)
      ag->primary[iPr].c[i] = tmp[i] + ag->com[i];
  }
}

/******************************************************************************/
