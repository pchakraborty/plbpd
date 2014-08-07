#include "output.hpp"
#include "comm.hpp"
#include "error.hpp"
#include "lattice.hpp"
#include "domain.hpp"
#include "boundary.hpp"
#include "model.hpp"
#include "lb_dynamics.hpp"
#include "th_dynamics.hpp"
#include "p_dynamics.hpp"
#include "timer.hpp"
#include "input.hpp"
#include "restart.hpp"

#include "mpi.h"
#include <iomanip>
using std::ios;
using std::setw;
#include <cmath> // for fabs
#include <hdf5.h>

/******************************************************************************/

Output::Output(PLBPD* plbpd): Pointers(plbpd){
  forceFile = "output.F_hyd";
  torqueFile = "output.T_hyd";
  prmryTrajFile = "output.traj.xyz";
  aglmrtVelFile = "output.aggvel";
}

/******************************************************************************/

Output::~Output(){
  if(comm->root){
    if(Fout) Fout.close();
    if(Tout) Tout.close();
    if(tout) tout.close();
    if(vout) vout.close();
    if(logfile) fclose(logfile);
  }
}

/******************************************************************************/

void Output::setup(){

  // open force and torque output files and agglomerate trajectory file
  if(p_dynamics->nAglmrts>0){
    if(comm->root){
      // delete contents of the file if present
      Fout.open(forceFile.c_str(), std::ios::trunc);
      Fout.close();
      Fout.open(forceFile.c_str(), std::ios::app);
      if(!Fout){
	string str("could not open output file " + forceFile);
	error->all(str);
      }
      Tout.open(torqueFile.c_str(), std::ios::trunc);
      Tout.close();
      Tout.open(torqueFile.c_str(), std::ios::app);
      if(!Tout){
	string str("could not open output file " + torqueFile);
	error->all(str);
      }
      if(!input->aglmrt_fixed){
	tout.open(prmryTrajFile.c_str(), std::ios::trunc);
	tout.close();
	tout.open(prmryTrajFile.c_str(), std::ios::app);
	if(!tout){
	  string str("could not open output file " + prmryTrajFile);
	  error->all(str);
	}
	vout.open(aglmrtVelFile.c_str(), std::ios::trunc);
	vout.close();
	vout.open(aglmrtVelFile.c_str(), std::ios::app);
	if(!vout){
	  string str("could not open output file " + aglmrtVelFile);
	  error->all(str);
	}
      }
    }
  }

  // open log file
  if(comm->root){
    logfile = fopen("output.log","w");
    if(logfile==NULL){
      string str("could not open log file output.log");
      error->all(str);
    }
  }
}

/******************************************************************************/
  
void Output::write_moments(int tStep){
  MPI_Barrier(world);

  // create file in parallel
  char outFile[100];
  sprintf(outFile,"%s%d","output.moments.",tStep);
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,world,comm->info);
  hid_t file_id = H5Fcreate(outFile,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
  H5Pclose(plist_id);

  // call the approproate function
  restart->write_moments(file_id,plist_id);
  
  // close file
  H5Fclose(file_id);
}

/******************************************************************************/

void Output::write_force(int tStep){
  if(comm->root){
    for(unsigned int iAg=0; iAg<p_dynamics->nAglmrts; iAg++)
      Fout<<iAg<<" "<<tStep<<" "
	  <<p_dynamics->aglmrt[iAg].F_hyd[0]<<" "
	  <<p_dynamics->aglmrt[iAg].F_hyd[1]<<" "
	  <<p_dynamics->aglmrt[iAg].F_hyd[2]<<"\n";
    Fout.flush();
  }
}

/******************************************************************************/

void Output::write_torque(int tStep){
  if(comm->root){
    for(unsigned int iAg=0; iAg<p_dynamics->nAglmrts; iAg++)
      Tout<<iAg<<" "<<tStep<<" "
	  <<p_dynamics->aglmrt[iAg].T_hyd[0]<<" "
	  <<p_dynamics->aglmrt[iAg].T_hyd[1]<<" "
	  <<p_dynamics->aglmrt[iAg].T_hyd[2]<<"\n";
    Tout.flush();
  }
}

/******************************************************************************/

void Output::header(){
  MPI_Barrier(world);
  if(comm->root){
    fprintf(stdout,"\n");
    fprintf(logfile,"\n");
    if(input->lb_model=="bgk"){
      fprintf(stdout,"             LB model: BGK\n\n");
      fprintf(logfile,"             LB model: BGK\n\n");
    }
    else if(input->lb_model=="mrt"){
      if(input->fluctuations){
	fprintf(stdout,"             LB model: fluctuating MRT\n\n");
	fprintf(logfile,"             LB model: fluctuating MRT\n\n");
      }
      else{
	fprintf(stdout,"             LB model: MRT\n\n");
	fprintf(logfile,"             LB model: MRT\n\n");
      }
    }
    if(input->integrator=="rk4"){
      fprintf(stdout," integrator: rk4\n\n");
      fprintf(logfile," integrator: rk4\n\n");
    }
    else if(input->integrator=="rkf45"){
      fprintf(stdout," integrator: rkf45\n\n");
      fprintf(logfile," integrator: rkf45\n\n");
    }

    fprintf(stdout,"    time      avg(fluid rho)");
    if(input->heat_xfer) fprintf(stdout,"        avg(theta)\n");
    else fprintf(stdout,"\n");
    fprintf(logfile,"    time      avg(fluid rho)");
    if(input->heat_xfer) fprintf(logfile,"        avg(theta)\n");
    else fprintf(logfile,"\n");
    fflush(stdout);
    fflush(logfile);
  }
}

/******************************************************************************/

void Output::print2screen(int tStep){
  MPI_Barrier(world);
  double fluidrhoavg = lb_dynamics->calc_fluidrhoavg();
  double avgtemperature = 0.;
  if(input->heat_xfer)
    avgtemperature = th_dynamics->calc_avgtemperature();
  if(fabs(fluidrhoavg)>10.){
    string str("oops! abs[avg(rho)] for fluid is > 10.0");
    error->all(str);
  }
  if(isnan(fluidrhoavg)){
    string str = "oops! abs[avg(rho)] for system is 'nan'";
    error->all(str);
  }
  if(fabs(avgtemperature)>10.){
    string str = "oops! abs[avg(normalized temperature)] for system is > 10.0";
    error->all(str);
  }
  if(comm->root){
    fprintf(stdout,"%8i   %18.14f", tStep, fluidrhoavg);
    fprintf(logfile,"%8i   %18.14f", tStep, fluidrhoavg);
    if(input->heat_xfer){
      fprintf(stdout,"        %f", avgtemperature);
      fprintf(logfile,"        %f", avgtemperature);
    }
    fprintf(stdout,"\n");
    fprintf(logfile,"\n");
    
    fflush(stdout);
  }
}

/******************************************************************************/

void Output::write_nodetype_xz(){

  MPI_Barrier(world);
  int ***nodetype = lattice->type;

  char outFile[30];
  sprintf(outFile,"bdryxz.%d.dat",comm->me);
  boutxz.open(outFile);
    
  int yl = domain->Ly/2;
  for(int zl=0; zl<domain->lz+4; zl++){
    for(int xl=0; xl<domain->lx+4; xl++)
      boutxz<<setw(2)<<nodetype[zl][yl][xl]<<" ";
    boutxz<<"\n";
  }
  boutxz.close();
}

/******************************************************************************/

void Output::write_prmry_xyz(int tStep){
  MPI_Barrier(world);
  if(comm->root){
    double *cP=NULL;

    // total number of primaries
    int nP = 0;
    for(unsigned int iAg=0; iAg<p_dynamics->nAglmrts; iAg++)
      nP += p_dynamics->aglmrt[iAg].primary.size();

    tout<<" "<<nP<<"\n";
    tout<<" "<<tStep<<"\n";

    for(unsigned int iAg=0; iAg<p_dynamics->nAglmrts; iAg++){
      for(unsigned int iPr=0;iPr<p_dynamics->aglmrt[iAg].primary.size();iPr++){
	cP = p_dynamics->aglmrt[iAg].primary[iPr].c;
	tout<<" Al"<<" "<<cP[0]<<" "<<cP[1]<<" "<<cP[2]<<"\n";
      }
    }
    tout.flush();
  }
}

/******************************************************************************/

void Output::write_aglmrt_vel(int tStep){
  MPI_Barrier(world);
  if(comm->root){
    double *u = NULL;
    for(unsigned int iAg=0; iAg<p_dynamics->nAglmrts; iAg++){
      u = p_dynamics->aglmrt[iAg].u;
      vout<<tStep<<" "<<u[0]<<" "<<u[1]<<" "<<u[2]<<"\n";
    }
    vout.flush();
  }
}

/******************************************************************************/

void Output::nodecount(){
  
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  int intfld = 0;
  int inflow = 0;
  int outflow = 0;
  int noslip = 0;
  int freeslip = 0;
  int aglmrt = 0;

  int ***nodetype = lattice->type;

  // count only the domain nodes
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	if(nodetype[zl][yl][xl]==INTFLD)
	  intfld++;
	else if(nodetype[zl][yl][xl]==INFLOW)
	  inflow++;
	else if(nodetype[zl][yl][xl]==OUTFLOW)
	  outflow++;
	else if(nodetype[zl][yl][xl]==NOSLIP)
	  noslip++;
	else if(nodetype[zl][yl][xl]==FREESLIP)
	  freeslip++;
	else if(nodetype[zl][yl][xl]>=0)
	  aglmrt++;
	else{
	  string str("node type not recognized");
	  error->all(str);
	}
      }
    }
  }
  
  // add nodes across all procs
  int tmp;

  tmp = intfld;
  MPI_Reduce(&tmp,&intfld,1,MPI_INT,MPI_SUM,0,world);

  tmp = inflow;
  MPI_Reduce(&tmp,&inflow,1,MPI_INT,MPI_SUM,0,world);

  tmp = outflow;
  MPI_Reduce(&tmp,&outflow,1,MPI_INT,MPI_SUM,0,world);

  tmp = noslip;
  MPI_Reduce(&tmp,&noslip,1,MPI_INT,MPI_SUM,0,world);

  tmp = freeslip;
  MPI_Reduce(&tmp,&freeslip,1,MPI_INT,MPI_SUM,0,world);

  tmp = aglmrt;
  MPI_Reduce(&tmp,&aglmrt,1,MPI_INT,MPI_SUM,0,world);
  
  MPI_Barrier(world);
  if(comm->root){
    fprintf(stdout,"\n");
    fprintf(stdout," interior fluid nodes: %7d\n", intfld);
    fprintf(logfile," interior fluid nodes: %7d\n", intfld);
    fprintf(stdout,"   inflow fluid nodes: %7d\n", inflow);
    fprintf(logfile,"   inflow fluid nodes: %7d\n", inflow);
    fprintf(stdout,"  outflow fluid nodes: %7d\n", outflow);
    fprintf(logfile,"  outflow fluid nodes: %7d\n", outflow);
    fprintf(stdout,"   noslip solid nodes: %7d\n", noslip);
    fprintf(logfile,"   noslip solid nodes: %7d\n", noslip);
    fprintf(stdout," freeslip solid nodes: %7d\n", freeslip);
    fprintf(logfile," freeslip solid nodes: %7d\n", freeslip);
    fprintf(stdout,"   solid aglmrt nodes: %7d\n", aglmrt);
    fprintf(logfile,"   solid aglmrt nodes: %7d\n", aglmrt);
    tmp = intfld+inflow+outflow+noslip+freeslip+aglmrt;
    fprintf(stdout,"          total nodes: %7d\n", tmp);
    fprintf(logfile,"          total nodes: %7d\n", tmp);
    
    // a small check
    int nNodes = domain->Lx*domain->Ly*domain->Lz;
    if(nNodes!=tmp){
      string str("num of nodes don't match");
      error->all(str);
    }
  } 

}

/******************************************************************************/
void Output::write_time_taken(){
  MPI_Barrier(world);

  if(comm->root){
    double ml = timer->elapsed[MAINLOOP];
    double cNs_fl = timer->elapsed[COLLIDENSTREAM_FL];
    double abc_fl = timer->elapsed[APPLYBC_FL];
    double cm_fl = timer->elapsed[CALCMOMENTS_FL];
    double cNs_th = timer->elapsed[COLLIDENSTREAM_TH];
    double abc_th = timer->elapsed[APPLYBC_TH];
    double cm_th = timer->elapsed[CALCMOMENTS_TH];
    double comm = timer->elapsed[COMM];
    double out = timer->elapsed[OUTPUT];
    double mdi = timer->elapsed[INTEGRATE];
    double rfm = timer->elapsed[RESET_FLAG_MOMENT];

    fprintf(stdout,"\n times taken...\n\n");
    fprintf(logfile,"\n times taken...\n\n");

    fprintf(stdout,"   main loop: %f\n", ml);
    fprintf(logfile,"   main loop: %f\n", ml);

    fprintf(stdout,"   fluid dynamics:\n");
    fprintf(logfile,"   fluid dynamics:\n");
    fprintf(stdout,"     collide and stream: %f (%05.2f%%)\n",
	    cNs_fl, cNs_fl*100.0/ml);
    fprintf(logfile,"     collide and stream: %f (%05.2f%%)\n",
	    cNs_fl, cNs_fl*100.0/ml);

    fprintf(stdout,"     apply bdry condn: %f (%05.2f%%)\n",
	    abc_fl, abc_fl*100.0/ml);
    fprintf(logfile,"     apply bdry condn: %f (%05.2f%%)\n",
	    abc_fl, abc_fl*100.0/ml);

    fprintf(stdout,"     calculate moments: %f (%05.2f%%)\n",
	    cm_fl, cm_fl*100.0/ml);
    fprintf(logfile,"     calculate moments: %f (%05.2f%%)\n",
	    cm_fl, cm_fl*100.0/ml);

    if(input->heat_xfer){
      fprintf(stdout,"   heat transfer:\n");
      fprintf(logfile,"   heat transfer:\n");
      fprintf(stdout,"     collide and stream: %f (%05.2f%%)\n",
	      cNs_th, cNs_th*100.0/ml);
      fprintf(logfile,"     collide and stream: %f (%05.2f%%)\n",
	      cNs_th, cNs_th*100.0/ml);
      
      fprintf(stdout,"     apply bdry condn: %f (%05.2f%%)\n",
	      abc_th, abc_th*100.0/ml);
      fprintf(logfile,"     apply bdry condn: %f (%05.2f%%)\n",
	      abc_th, abc_th*100.0/ml);
      
      fprintf(stdout,"     calculate moments: %f (%05.2f%%)\n",
	      cm_th, cm_th*100.0/ml);
      fprintf(logfile,"     calculate moments: %f (%05.2f%%)\n",
	      cm_th, cm_th*100.0/ml);
    }

    if(input->nAg>0){
      fprintf(stdout,"   particle dynamics:\n");
      fprintf(logfile,"   particle dynamics:\n");
      fprintf(stdout,"     pd integrate: %f (%05.2f%%)\n", mdi, mdi*100.0/ml);
      fprintf(logfile,"     pd integrate: %f (%05.2f%%)\n", mdi, mdi*100.0/ml);
      fprintf(stdout,"     reset flag: %f (%05.2f%%)\n", rfm, rfm*100.0/ml);
      fprintf(logfile,"     reset flag: %f (%05.2f%%)\n", rfm, rfm*100.0/ml);
    }
    
    fprintf(stdout,"   communication: %f (%05.2f%%)\n", comm, comm*100.0/ml);
    fprintf(logfile,"   communication: %f (%05.2f%%)\n", comm, comm*100.0/ml);

    fprintf(stdout,"   output: %f (%05.2f%%)\n", out, out*100.0/ml);
    fprintf(logfile,"   output: %f (%05.2f%%)\n", out, out*100.0/ml);

    fprintf(stdout,"\n ...done!\n");
    fprintf(logfile,"\n ...done!\n");
  }  
}

/******************************************************************************/

void Output::write_restart(int tStep){
  MPI_Barrier(world);
  
  // create file in parallel
  char outFile[100];
  sprintf(outFile,"%s%d","output.restart.",tStep);
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,world,comm->info);
  hid_t file_id = H5Fcreate(outFile,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
  H5Pclose(plist_id);

  restart->write_restart(file_id, plist_id, tStep);

  // close file
  H5Fclose(file_id);
}

/******************************************************************************/
