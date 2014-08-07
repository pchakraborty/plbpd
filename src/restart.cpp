#include "restart.hpp"
#include "input.hpp"
#include "lattice.hpp"
#include "model.hpp"
#include "domain.hpp"
#include "comm.hpp"
#include "p_dynamics.hpp"
#include "agglomerate.hpp"
#include "memory.hpp"
#include "primary.hpp"
#include "error.hpp"
#include <mpi.h>
#include <cstring>

/******************************************************************************/

Restart::Restart(PLBPD *plbpd): Pointers(plbpd){
  inp_params_id = 0;
  filespace_inp = 0;
  agdata_id = 0;

  filespace = 0;
  memspace = 0;

  filespace_3 = 0;
  memspace_3 = 0;

  filespace_nV_fl = 0;
  memspace_nV_fl = 0;

  filespace_nV_th = 0;
  memspace_nV_th = 0;
}

/******************************************************************************/

Restart::~Restart(){
  if(inp_params_id>0) H5Tclose(inp_params_id);
  if(filespace_inp>0) H5Sclose(filespace_inp);
  if(agdata_id>0) H5Tclose(agdata_id);

  if(filespace>0) H5Sclose(filespace);
  if(memspace>0) H5Sclose(memspace);

  if(filespace_3>0) H5Sclose(filespace_3);
  if(memspace_3>0) H5Sclose(memspace_3);

  if(filespace_nV_fl>0) H5Sclose(filespace_nV_fl);
  if(memspace_nV_fl>0) H5Sclose(memspace_nV_fl);

  if(filespace_nV_th>0) H5Sclose(filespace_nV_th);
  if(memspace_nV_th>0) H5Sclose(memspace_nV_th);
}

/******************************************************************************/

void Restart::setup(){
  // input variables
  {
    hid_t rank = 1;
    hsize_t dim[] = {1};
    filespace_inp = H5Screate_simple(rank,dim,NULL);

    // create compound datatype for memory
    hsize_t ARR3_SIZE[] = {3};
    hid_t ARR3_DBL = H5Tarray_create1(H5T_NATIVE_DOUBLE,1,ARR3_SIZE,NULL);
    hsize_t ARR6_SIZE[] = {6};
    hid_t ARR6_CHAR = H5Tarray_create1(H5T_NATIVE_CHAR,1,ARR6_SIZE,NULL);
    hid_t ARR6_DBL = H5Tarray_create1(H5T_NATIVE_DOUBLE,1,ARR6_SIZE,NULL);
    hsize_t ARR10_SIZE[] = {10};
    hid_t ARR10_CHAR = H5Tarray_create1(H5T_NATIVE_CHAR,1,ARR10_SIZE,NULL);

    inp_params_id = H5Tcreate(H5T_COMPOUND, sizeof(inp_params));
    H5Tinsert(inp_params_id,"tstep",HOFFSET(inp_params,tstep),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"npx",HOFFSET(inp_params,npx),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"npy",HOFFSET(inp_params,npy),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"npz",HOFFSET(inp_params,npz),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"Lx",HOFFSET(inp_params,Lx),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"Ly",HOFFSET(inp_params,Ly),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"Lz",HOFFSET(inp_params,Lz),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"bdryType",HOFFSET(inp_params,bdryType),ARR6_CHAR);
    H5Tinsert(inp_params_id,"lb_model",HOFFSET(inp_params,lb_model),ARR10_CHAR);
    H5Tinsert(inp_params_id,"nu",HOFFSET(inp_params,nu),H5T_NATIVE_DOUBLE);
    H5Tinsert(inp_params_id,"uE",HOFFSET(inp_params,uE),ARR3_DBL);
    H5Tinsert(inp_params_id,"uW",HOFFSET(inp_params,uW),ARR3_DBL);
    H5Tinsert(inp_params_id,"uN",HOFFSET(inp_params,uN),ARR3_DBL);
    H5Tinsert(inp_params_id,"uS",HOFFSET(inp_params,uS),ARR3_DBL);
    H5Tinsert(inp_params_id,"uU",HOFFSET(inp_params,uU),ARR3_DBL);
    H5Tinsert(inp_params_id,"uD",HOFFSET(inp_params,uD),ARR3_DBL);
    H5Tinsert(inp_params_id,"extF",HOFFSET(inp_params,extF),ARR3_DBL);
    // for mrt
    H5Tinsert(inp_params_id,"nu_b",HOFFSET(inp_params,nu_b),H5T_NATIVE_DOUBLE);
    H5Tinsert(inp_params_id,"temperature",
	      HOFFSET(inp_params,temperature),H5T_NATIVE_DOUBLE);
    H5Tinsert(inp_params_id,"fluct_seed",
	      HOFFSET(inp_params,fluct_seed),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"nAg",HOFFSET(inp_params,nAg),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"aglmrt_fixed",
	      HOFFSET(inp_params,aglmrt_fixed),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"p_start_step",
	      HOFFSET(inp_params,p_start_step),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"integrator",
	      HOFFSET(inp_params,integrator),ARR10_CHAR);
    // for non-isothermal flow
    H5Tinsert(inp_params_id,"heat_xfer",
	      HOFFSET(inp_params,heat_xfer),H5T_NATIVE_INT);
    H5Tinsert(inp_params_id,"Ra",
	      HOFFSET(inp_params,Ra),H5T_NATIVE_DOUBLE);
    H5Tinsert(inp_params_id,"th_diff_fl",
	      HOFFSET(inp_params,th_diff_fl),H5T_NATIVE_DOUBLE);
    H5Tinsert(inp_params_id,"bdryType_th",
	      HOFFSET(inp_params,bdryType_th),ARR6_CHAR); 
    H5Tinsert(inp_params_id,"bdryTemp",
	      HOFFSET(inp_params,bdryTemp),ARR6_DBL); 
  }
  
  // agglomerate data
  if(input->nAg>0){
    hid_t rank = 1;
    hsize_t dim[] ={input->nAg};
    filespace_ag = H5Screate_simple(rank,dim,NULL);

    // create compound datatype for memory
    hsize_t ARR3_SIZE[] = {3};
    hid_t ARR3_TYPE = H5Tarray_create1(H5T_NATIVE_DOUBLE,1,ARR3_SIZE,NULL);
    hsize_t ARR9_SIZE[] = {9};
    hid_t ARR9_TYPE = H5Tarray_create1(H5T_NATIVE_DOUBLE,1,ARR9_SIZE,NULL);
    hsize_t ARRN_SIZE[] = {PLBPD_NS::MAXPRMRYS*PLBPD_NS::NPRMRYDBLES};
    hid_t ARRN_TYPE = H5Tarray_create1(H5T_NATIVE_DOUBLE,1,ARRN_SIZE,NULL);

    agdata_id = H5Tcreate(H5T_COMPOUND, sizeof(ag_data));
    H5Tinsert(agdata_id,"ID",      HOFFSET(ag_data,ID),H5T_NATIVE_UINT);
    H5Tinsert(agdata_id,"nPr",     HOFFSET(ag_data,nPr),H5T_NATIVE_UINT);
    H5Tinsert(agdata_id,"prmry",   HOFFSET(ag_data,prmry),ARRN_TYPE);
    H5Tinsert(agdata_id,"rho",     HOFFSET(ag_data,rho),H5T_NATIVE_DOUBLE);
    H5Tinsert(agdata_id,"vol",     HOFFSET(ag_data,vol),H5T_NATIVE_DOUBLE);
    H5Tinsert(agdata_id,"th_diff", HOFFSET(ag_data,th_diff),H5T_NATIVE_DOUBLE);
    H5Tinsert(agdata_id,"mass",    HOFFSET(ag_data,mass),H5T_NATIVE_DOUBLE);
    H5Tinsert(agdata_id,"Ibody",   HOFFSET(ag_data,Ib),ARR3_TYPE);
    H5Tinsert(agdata_id,"IbodyInv",HOFFSET(ag_data,IbI),ARR3_TYPE);
    H5Tinsert(agdata_id,"com",     HOFFSET(ag_data,com),ARR3_TYPE);
    H5Tinsert(agdata_id,"RotMat",  HOFFSET(ag_data,R),ARR9_TYPE);
    H5Tinsert(agdata_id,"LinMom",  HOFFSET(ag_data,P),ARR3_TYPE);
    H5Tinsert(agdata_id,"AngMom",  HOFFSET(ag_data,L),ARR3_TYPE);
    H5Tinsert(agdata_id,"Iinv",    HOFFSET(ag_data,Iinv),ARR9_TYPE);
    H5Tinsert(agdata_id,"u",       HOFFSET(ag_data,u),ARR3_TYPE);
    H5Tinsert(agdata_id,"omega",   HOFFSET(ag_data,omega),ARR3_TYPE);
    H5Tinsert(agdata_id,"F_hyd",   HOFFSET(ag_data,F_hyd),ARR3_TYPE);
    H5Tinsert(agdata_id,"T_hyd",   HOFFSET(ag_data,T_hyd),ARR3_TYPE);
    H5Tinsert(agdata_id,"F_tot",   HOFFSET(ag_data,F_tot),ARR3_TYPE);
    H5Tinsert(agdata_id,"T_tot",   HOFFSET(ag_data,T_tot),ARR3_TYPE);
  }
}

/******************************************************************************/

void Restart::alloc(){
  // allocate memory for data & memory spaces for binary output
  int Lz = domain->Lz;
  int Ly = domain->Ly;
  int Lx = domain->Lx;
  int lz = domain->lz;
  int ly = domain->ly;
  int lx = domain->lx;
  int lzmin = domain->lzmin;
  int lymin = domain->lymin;
  int lxmin = domain->lxmin;

  { // for type/rho/alpha/theta
    hid_t rank = 3;
    hsize_t dim[3], offset[3];
    hsize_t count[] = {lz,ly,lx};
    // filespace
    dim[0] = Lz; dim[1] = Ly; dim[2] = Lx;
    filespace = H5Screate_simple(rank,dim,NULL);
    offset[0] = lzmin; offset[1] = lymin; offset[2] = lxmin;
    H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);
    // memoryspace
    dim[0] = lz+2*2; dim[1] = ly+2*2; dim[2] = lx+2*2;
    memspace = H5Screate_simple(rank,dim,NULL);
    offset[0] = 2; offset[1] = 2; offset[2] = 2;
    H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,NULL,count,NULL);
  }
    
  { // for u[3]
    hid_t rank = 4;
    hsize_t dim[4], offset[4];
    hsize_t count[] = {lz,ly,lx,3};
    // filespace
    dim[0] = Lz; dim[1] = Ly; dim[2] = Lx; dim[3] = 3;
    filespace_3 = H5Screate_simple(rank,dim,NULL);
    offset[0] = lzmin; offset[1] = lymin; offset[2] = lxmin; offset[3] = 0;
    H5Sselect_hyperslab(filespace_3,H5S_SELECT_SET,offset,NULL,count,NULL);
    // memoryspace
    dim[0] = lz+2*2; dim[1] = ly+2*2; dim[2] = lx+2*2; dim[3] = 3;
    memspace_3 = H5Screate_simple(rank,dim,NULL);
    offset[0] = 2; offset[1] = 2; offset[2] = 2; offset[3] = 0;
    H5Sselect_hyperslab(memspace_3,H5S_SELECT_SET,offset,NULL,count,NULL);
  }    

  { // for n[nVel] of fluid
    hid_t rank = 4;
    hsize_t dim[4], offset[4];
    int nVel_fl = model->fluid->nVelocity;
    hsize_t count[] = {lz,ly,lx,nVel_fl};
    // filespace
    dim[0] = Lz; dim[1] = Ly; dim[2] = Lx; dim[3] = nVel_fl;
    filespace_nV_fl = H5Screate_simple(rank,dim,NULL);
    offset[0] = lzmin; offset[1] = lymin; offset[2] = lxmin; offset[3] = 0;
    H5Sselect_hyperslab(filespace_nV_fl,H5S_SELECT_SET,offset,NULL,count,NULL);
    // memoryspace
    dim[0] = lz+2*2; dim[1] = ly+2*2; dim[2] = lx+2*2; dim[3] = nVel_fl;
    memspace_nV_fl = H5Screate_simple(rank,dim,NULL);
    offset[0] = 2; offset[1] = 2; offset[2] = 2; offset[3] = 0;
    H5Sselect_hyperslab(memspace_nV_fl,H5S_SELECT_SET,offset,NULL,count,NULL);
  }  
  
  // for n[nVel] of heat
  if(input->heat_xfer){
    hid_t rank = 4;
    hsize_t dim[4], offset[4];
    int nVel_th = model->thermal->nVelocity;
    hsize_t count[] = {lz,ly,lx,nVel_th};
    // filespace
    dim[0] = Lz; dim[1] = Ly; dim[2] = Lx; dim[3] = nVel_th;
    filespace_nV_th = H5Screate_simple(rank,dim,NULL);
    offset[0] = lzmin; offset[1] = lymin; offset[2] = lxmin; offset[3] = 0;
    H5Sselect_hyperslab(filespace_nV_th,H5S_SELECT_SET,offset,NULL,count,NULL);
    // memoryspace
    dim[0] = lz+3*2; dim[1] = ly+3*2; dim[2] = lx+3*2; dim[3] = nVel_th;
    memspace_nV_th = H5Screate_simple(rank,dim,NULL);
    offset[0] = 3; offset[1] = 3; offset[2] = 3; offset[3] = 0;
    H5Sselect_hyperslab(memspace_nV_th,H5S_SELECT_SET,offset,NULL,count,NULL);
  }      
}

/******************************************************************************/

void Restart::write_moments(hid_t file_id, hid_t plist_id){
  
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE); // COLLECTIVE
  
  // write rho
  {
    // create dataset
    hid_t dataset;
    dataset = H5Dcreate1(file_id,"rho",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT);
    
    // write to disk
    double *ptr = &lattice->rho[0][0][0];
    H5Dwrite(dataset,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,ptr);
    
    // release resources
    H5Dclose(dataset);
  }
  // write u
  {
    // create dataset
    hid_t dataset;
    dataset = H5Dcreate1(file_id,"u",H5T_NATIVE_DOUBLE,filespace_3,H5P_DEFAULT);
    
    // write to disk
    double *ptr = &lattice->u[0][0][0][0];
    H5Dwrite(dataset,H5T_NATIVE_DOUBLE,memspace_3,filespace_3,plist_id,ptr);
    
    // release resources
    H5Dclose(dataset);
  }
  // write temperature (in case of non-isothermal flow)
  if(input->heat_xfer){
    // create dataset
    hid_t dataset;
    dataset=H5Dcreate1(file_id,"theta",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT);
    
    // write to disk
    double *ptr = &lattice->theta[0][0][0];
    H5Dwrite(dataset,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,ptr);
    
    // release resources
    H5Dclose(dataset);
  }
  
  // release resources
  H5Pclose(plist_id);
}

/******************************************************************************/

void Restart::write_restart(hid_t file_id, hid_t plist_id, int timestep){

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE); // COLLECTIVE

  // write important input parameters
  { // only root writes input parameters
    inp_params i;

    // some initializations, so that we don't have
    // numbers like 1.52063e-314 in the restart file
    for(int j=0; j<10; j++){
      i.lb_model[j] = 0;
      i.integrator[j] = 0;
    }
    for(int j=0; j<6; j++){
      i.bdryTemp[j] = 0.0;
      i.bdryType_th[j] = 0;
    }
    i.Ra = -1.0;
    i.th_diff_fl = -1.0;
    i.nu_b = 0.;
    i.temperature = 0.;
    i.fluct_seed = 0;
    i.aglmrt_fixed = 1;
    i.p_start_step = 0;

    i.tstep = timestep;
    i.npx = input->npx; i.npy = input->npy; i.npz = input->npz;
    i.Lx = input->Lx; i.Ly = input->Ly; i.Lz = input->Lz;
    for(int j=0; j<6; j++)
      i.bdryType[j] = input->bdryType[j][0];
    strcpy(i.lb_model,(input->lb_model).c_str());
    i.nu = input->nu;
    for(int j=0; j<3; j++){
      i.uE[j] = input->uE[j];
      i.uW[j] = input->uW[j];
      i.uN[j] = input->uN[j];
      i.uS[j] = input->uS[j];
      i.uU[j] = input->uU[j];
      i.uD[j] = input->uD[j];
      i.extF[j] = input->extF[j];
    }
    // mrt
    if(input->lb_model=="mrt"){
      i.nu_b = input->nu_b;
      i.temperature = input->temperature;
      i.fluct_seed = input->fluct_seed;
    }
    // aglmrts
    i.nAg = input->nAg;
    if(i.nAg>0){
      if(input->aglmrt_fixed)
	i.aglmrt_fixed = 1; // True
      else
	i.aglmrt_fixed = 0; // False
      i.p_start_step = input->p_start_step;
      strcpy(i.integrator,(input->integrator).c_str());
    }
    // heat transfer
    if(input->heat_xfer){
      i.heat_xfer = 1;    // True
      i.Ra = input->Ra;
      i.th_diff_fl = input->th_diff_fl;
      for(int j=0; j<6; j++){
	i.bdryType_th[j] = input->bdryType_th[j][0];
	i.bdryTemp[j] = input->bdryTemp[j];
      }
    }
    else{
      i.heat_xfer = 0;    // False
    }
    
    // create dataset
    hid_t dataset;
    dataset = H5Dcreate1(file_id,"params",inp_params_id,
			filespace_inp,H5P_DEFAULT);
    
    // write to file
    if(comm->root)
      H5Dwrite(dataset,inp_params_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,&i);
    
    // release resources
    H5Dclose(dataset);
  }

  // only root writes agglomerate data
  if(input->nAg>0){
    Agglomerate *aglmrt = p_dynamics->aglmrt;
    ag_data agdata[input->nAg];

    for(int iAg=0; iAg<input->nAg; iAg++){
      agdata[iAg].ID = aglmrt[iAg].ID;
      agdata[iAg].nPr = aglmrt[iAg].nPrimaries;
      for(unsigned int i=0; i<PLBPD_NS::MAXPRMRYS*PLBPD_NS::NPRMRYDBLES; i++)
	agdata[iAg].prmry[i] = 0.;
      for(unsigned int iPr=0; iPr<aglmrt[iAg].nPrimaries; iPr++){
	agdata[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+0] = 
	  aglmrt[iAg].primary[iPr].r;
	agdata[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+1] = 
	  aglmrt[iAg].primary[iPr].c[0];
	agdata[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+2] = 
	  aglmrt[iAg].primary[iPr].c[1];
	agdata[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+3] =
	  aglmrt[iAg].primary[iPr].c[2];
	agdata[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+4] =
	  aglmrt[iAg].primary[iPr].cb[0];
	agdata[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+5] =
	  aglmrt[iAg].primary[iPr].cb[1];
	agdata[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+6] =
	  aglmrt[iAg].primary[iPr].cb[2];
      }
      agdata[iAg].rho = aglmrt[iAg].rho;
      agdata[iAg].vol = aglmrt[iAg].volume;
      if(input->heat_xfer)
	agdata[iAg].th_diff = aglmrt[iAg].th_diff;
      agdata[iAg].mass = aglmrt[iAg].mass;
      for(int i=0; i<3; i++){
	agdata[iAg].Ib[i] = aglmrt[iAg].Ibody[i];
	agdata[iAg].IbI[i] = aglmrt[iAg].IbodyInv[i];
	agdata[iAg].com[i] = aglmrt[iAg].com[i];
	agdata[iAg].P[i] = aglmrt[iAg].P[i];
	agdata[iAg].L[i] = aglmrt[iAg].L[i];
	agdata[iAg].u[i] = aglmrt[iAg].u[i];
	agdata[iAg].omega[i] = aglmrt[iAg].omega[i];
	agdata[iAg].F_hyd[i] = aglmrt[iAg].F_hyd[i];
	agdata[iAg].T_hyd[i] = aglmrt[iAg].T_hyd[i];
	agdata[iAg].F_tot[i] = aglmrt[iAg].F_tot[i];
	agdata[iAg].T_tot[i] = aglmrt[iAg].T_tot[i];
      }
      for(int i=0; i<9; i++){
	agdata[iAg].R[i] = aglmrt[iAg].R[i];
	agdata[iAg].Iinv[i] = aglmrt[iAg].Iinv[i];
      }
    }
    
    // create dataset
    hid_t dataset;
    dataset = H5Dcreate1(file_id,"aglmrt",agdata_id,filespace_ag,H5P_DEFAULT);

    // write to file
    if(comm->root)
      H5Dwrite(dataset,agdata_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,agdata);

    // release resources
    H5Dclose(dataset);
  }


  // rest are written collectively
  // write n
  {
    // create dataset
    hid_t dataset;
    dataset=H5Dcreate1(file_id,"n",H5T_NATIVE_DOUBLE,
		      filespace_nV_fl,H5P_DEFAULT);

    // write to disk
    double *n = &lattice->n[0][0][0][0];
    H5Dwrite(dataset,H5T_NATIVE_DOUBLE,
	     memspace_nV_fl,filespace_nV_fl,plist_id,n);
    
    // release resources
    H5Dclose(dataset);
  }

  // write type
  {
    // create dataset
    hid_t dataset;
    dataset=H5Dcreate1(file_id,"type",H5T_NATIVE_INT,filespace,H5P_DEFAULT);
    
    // write to disk
    int *type = &lattice->type[0][0][0];
    H5Dwrite(dataset,H5T_NATIVE_INT,memspace,filespace,plist_id,type);

    // release resources
    H5Dclose(dataset);
  }

  // write moments
  write_moments(file_id, plist_id);

  // write temperature density for the case of non-isothermal flow
  if(input->heat_xfer){
    // create dataset
    hid_t dataset;
    dataset=H5Dcreate1(file_id,"nth",H5T_NATIVE_DOUBLE,
		      filespace_nV_th,H5P_DEFAULT);

    // write to disk
    double *nth = &lattice->nth[0][0][0][0];
    H5Dwrite(dataset,H5T_NATIVE_DOUBLE,
	     memspace_nV_th,filespace_nV_th,plist_id,nth);
    
    // release resources
    H5Dclose(dataset);
  }

  // write thermal diffusivity for the case of non-isothermal flow
  if(input->heat_xfer){
    // create dataset
    hid_t dataset;
    dataset=H5Dcreate1(file_id,"alpha",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT);
    
    // write to disk
    double *alpha = &lattice->alpha[0][0][0];
    H5Dwrite(dataset,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,alpha);
    
    // release resources
    H5Dclose(dataset);
  }

  // release resources
  H5Pclose(plist_id);
}

/******************************************************************************/

void Restart::read_inpparams(hid_t file_id){
  inp_params inpp;

  // open dataset 'params'
  hid_t dataset = H5Dopen1(file_id,"/params");

  // read input parameters
  H5Dread(dataset,inp_params_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,&inpp);

  // release resources
  H5Dclose(dataset);

  // copy the rest of the input parameters
  plbpd->timestep = inpp.tstep;
  input->Lx = inpp.Lx;
  input->Ly = inpp.Ly;
  input->Lz = inpp.Lz;
  input->nu = inpp.nu;
  for(int i=0; i<3; i++){
    input->uE[i] = inpp.uE[i];
    input->uW[i] = inpp.uW[i];
    input->uN[i] = inpp.uN[i];
    input->uS[i] = inpp.uS[i];
    input->uU[i] = inpp.uU[i];
    input->uD[i] = inpp.uD[i];
    input->extF[i] = inpp.extF[i];
  }
  // mrt
  input->nu_b = inpp.nu_b;
  input->temperature = inpp.temperature;
  // aglmrt
  if(inpp.aglmrt_fixed==1)
    input->aglmrt_fixed = true;
  else
    input->aglmrt_fixed = false;

  input->p_start_step = inpp.p_start_step;
  input->integrator = inpp.integrator;

  // for non-isothermal flpiows
  if(inpp.heat_xfer==1){
    input->heat_xfer = true;
    input->Ra = inpp.Ra;
    input->th_diff_fl = inpp.th_diff_fl;
    for(int i=0; i<6; i++)
      input->bdryTemp[i] = inpp.bdryTemp[i];
  }
  else if(inpp.heat_xfer==0)
    input->heat_xfer = false;
  else{
    string str = "input parameter 'heat_xfer' in restart file not understood";
    error->all(str);
  }

  // some values (e.g lb_model, Lx etc) are read both
  // from the input file and the restart input parameters.
  // check to see that if they are equal

  if((input->npx!=inpp.npx)or(input->npy!=inpp.npy)or(input->npz!=inpp.npz)){
    string str =
      "we are using a set of procs different from the original simulation";
    error->notify(str);
  }
  if(strcmp(input->lb_model.c_str(),inpp.lb_model)!=0){
    string str = "lb_model read from input file '" + input->lb_model +
      "'should be the same as the one read from the restart file" +
      inpp.lb_model;
    error->all(str);
  }
  if(input->nAg!=inpp.nAg){
    string str = "num_aglmrts read from the input file ("
      + utils::T2str<int>(input->nAg)
      + ") should be the same as the one read from the restart file (" 
      + utils::T2str<int>(inpp.nAg) + ")";
    error->all(str);
  }
  for(int j=0; j<6; j++){
    if(input->bdryType[j][0]!=inpp.bdryType[j]){
      string str = "bdryType read from restart file is ";
      str += "NOT the same as the ones read from the input file";
      error->all(str);
    }
    if(input->bdryType_th[j][0]!=inpp.bdryType_th[j]){
      string str = "bdryType_th read from restart file is NOT ";
      str += "the same as the ones read from the input file";
      error->all(str);
    }
  }

}

/******************************************************************************/

void Restart::read_data(hid_t file_id){
  // read rho
  {
    double ***rho = lattice->rho;

    // open dataset
    hid_t dataset = H5Dopen1(file_id,"rho");

    // read data
    H5Dread(dataset,H5T_NATIVE_DOUBLE,
	    memspace,filespace,H5P_DEFAULT,&rho[0][0][0]);

    // release resource
    H5Dclose(dataset);
  }

  // read type
  {
    int ***type = lattice->type;

    // open dataset
    hid_t dataset = H5Dopen1(file_id,"type");

    // read data
    H5Dread(dataset,H5T_NATIVE_INT,
	    memspace,filespace,H5P_DEFAULT,&type[0][0][0]);

    // release resource
    H5Dclose(dataset);
  }
    
  // read u
  {
    double ****u = lattice->u;

    // open dataset
    hid_t dataset = H5Dopen1(file_id,"u");

    // read data
    H5Dread(dataset,H5T_NATIVE_DOUBLE,
	    memspace_3,filespace_3,H5P_DEFAULT,&u[0][0][0][0]);

    // release resource
    H5Dclose(dataset);
  }

  // read n
  {
    double ****n = lattice->n;

    // open dataset
    hid_t dataset = H5Dopen1(file_id,"n");

    // read data
    H5Dread(dataset,H5T_NATIVE_DOUBLE,
	    memspace_nV_fl,filespace_nV_fl,H5P_DEFAULT,&n[0][0][0][0]);

    // release resource
    H5Dclose(dataset);
  }

  // read nth (for non-isothermal flow)
  if(input->heat_xfer){
    double ****nth = lattice->nth;

    // open dataset
    hid_t dataset = H5Dopen1(file_id,"nth");

    // read data
    H5Dread(dataset,H5T_NATIVE_DOUBLE,
	    memspace_nV_th,filespace_nV_th,H5P_DEFAULT,&nth[0][0][0][0]);

    // release resource
    H5Dclose(dataset);
  }
    
  // read theta (for non-isothermal flow)
  if(input->heat_xfer){
    double ***theta = lattice->theta;

    // open dataset
    hid_t dataset = H5Dopen1(file_id,"theta");

    // read data
    H5Dread(dataset,H5T_NATIVE_DOUBLE,
	    memspace,filespace,H5P_DEFAULT,&theta[0][0][0]);

    // release resources
    H5Dclose(dataset);
  }

  // read alpha (for non-isothermal flow)
  if(input->heat_xfer){
    double ***alpha = lattice->alpha;
    
    // open dataset
    hid_t dataset = H5Dopen1(file_id,"alpha");

    // read data
    H5Dread(dataset,H5T_NATIVE_DOUBLE,
	    memspace,filespace,H5P_DEFAULT,&alpha[0][0][0]);

    // release resources
    H5Dclose(dataset);
  }

  // finally read aglmrt data
  if(input->nAg>0){
    Agglomerate *aglmrt = p_dynamics->aglmrt;
    ag_data *ag = new ag_data[input->nAg];

    // open dataset
    hid_t dataset = H5Dopen1(file_id,"aglmrt");

    // read
    for(int iAg=0; iAg<input->nAg; iAg++)
      H5Dread(dataset,agdata_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,&ag[iAg]);

    // release resources
    H5Dclose(dataset);

    // copy
    int prmryID = 0;
    for(int iAg=0; iAg<input->nAg; iAg++){
      aglmrt[iAg].ID = ag[iAg].ID;
      aglmrt[iAg].nPrimaries = ag[iAg].nPr;
      aglmrt[iAg].rho = ag[iAg].rho;
      aglmrt[iAg].volume = ag[iAg].vol;
      if(input->heat_xfer)
	aglmrt[iAg].th_diff = ag[iAg].th_diff;
      aglmrt[iAg].mass = ag[iAg].mass;

      for(int j=0; j<3; j++){
	aglmrt[iAg].Ibody[j] = ag[iAg].Ib[j];
	aglmrt[iAg].IbodyInv[j] = ag[iAg].IbI[j];
	aglmrt[iAg].com[j] = ag[iAg].com[j];
	aglmrt[iAg].P[j] = ag[iAg].P[j];
	aglmrt[iAg].L[j] = ag[iAg].L[j];
	aglmrt[iAg].u[j] = ag[iAg].u[j];
	aglmrt[iAg].omega[j] = ag[iAg].omega[j];
	aglmrt[iAg].F_hyd[j] = ag[iAg].F_hyd[j];
	aglmrt[iAg].T_hyd[j] = ag[iAg].T_hyd[j];
	aglmrt[iAg].F_tot[j] = ag[iAg].F_tot[j];
	aglmrt[iAg].T_tot[j] = ag[iAg].T_tot[j];
      }
      for(int j=0; j<9; j++){
	aglmrt[iAg].R[j] = ag[iAg].R[j];
	aglmrt[iAg].Iinv[j] = ag[iAg].Iinv[j];
      }
      Primary tmp;
      for(unsigned int iPr=0; iPr<aglmrt[iAg].nPrimaries; iPr++){
	tmp.ID    = prmryID;
	tmp.r     = ag[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+0];
	tmp.c[0]  = ag[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+1];
	tmp.c[1]  = ag[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+2];
	tmp.c[2]  = ag[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+3];
	tmp.cb[0] = ag[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+4];
	tmp.cb[1] = ag[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+5];
	tmp.cb[2] = ag[iAg].prmry[iPr*PLBPD_NS::NPRMRYDBLES+6];

	aglmrt[iAg].primary.push_back(tmp);
	prmryID++;
      }
    }
  }
}

/******************************************************************************/
