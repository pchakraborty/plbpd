#include <hdf5.h>
#include <mpi.h>
#include "input.hpp"
#include "tuples.hpp"
#include "plbpd.hpp"
#include "restart.hpp"
#include "error.hpp"
#include "config.hpp"

/******************************************************************************/

Input::Input(PLBPD *plbpd): Pointers(plbpd){
}

/******************************************************************************/

Input::~Input(){
  fluctuations = false;
}

/******************************************************************************/

void Input::read(){
  /* open parameter file */
  config->open(plbpd->infile);
  
  tuples3<double> v(0.,0.,0.);

  // processors
  tuples3<int> procs(0,0,0);
  procs = config->read<tuples3<int> >("processors");
  npx = procs.a; npy = procs.b; npz = procs.c;

  // domain dimension
  if(!plbpd->RSTRT){
    tuples3<int> dim(0,0,0);
    dim = config->read<tuples3<int> >("domain_dim");
    Lx = dim.a; Ly = dim.b; Lz = dim.c;
  } 
  
  // bdry type
  tuples6<string> bt("","","","","","");
  bt = config->read<tuples6<string> >("bdry_type");
  bdryType[0] = bt.a; bdryType[1] = bt.b; bdryType[2] = bt.c;
  bdryType[3] = bt.d; bdryType[4] = bt.e; bdryType[5] = bt.f;

  // LB model
  lb_model = config->read<string>("lb_model");

  // heat transfer/thermal dynamics
  string answer = config->read<string>("heat_transfer");
  if(answer=="yes")
    heat_xfer = true;
  else if(answer=="no")
    heat_xfer = false;
  else{
    string str = "input param 'heat_transfer' can be one of yes/no";
    error->all(str);
  }
  if(heat_xfer){
    bt = config->read<tuples6<string> >("bdry_type_th");
    bdryType_th[0] = bt.a; bdryType_th[1] = bt.b; bdryType_th[2] = bt.c;
    bdryType_th[3] = bt.d; bdryType_th[4] = bt.e; bdryType_th[5] = bt.f;
  }

  if(!plbpd->RSTRT){
    // rho(fluid)
    rhoFluid = config->read<double>("rho_fluid");
    
    // rho(solid domain)
    rhoSD = config->read<double>("rho_solid_domain");
    
    // nu
    nu = config->read<double>("nu");

    // flow velocity
    v = config->read<tuples3<double> >("flow_velocity");
    flowVel[0] = v.a; flowVel[1] = v.b; flowVel[2] = v.c;

    // boundary velocity: E
    v = config->read<tuples3<double> >("E_bdry_vel");
    uE[0] = v.a; uE[1] = v.b; uE[2] = v.c;
    
    // boundary velocity: W
    v = config->read<tuples3<double> >("W_bdry_vel");
    uW[0] = v.a; uW[1] = v.b; uW[2] = v.c;
    
    // boundary velocity: N
    v = config->read<tuples3<double> >("N_bdry_vel");
    uN[0] = v.a; uN[1] = v.b; uN[2] = v.c;
    
    // boundary velocity: S
    v = config->read<tuples3<double> >("S_bdry_vel");
    uS[0] = v.a; uS[1] = v.b; uS[2] = v.c;
    
    // boundary velocity: U
    v = config->read<tuples3<double> >("U_bdry_vel");
    uU[0] = v.a; uU[1] = v.b; uU[2] = v.c;
    
    // boundary velocity: D
    v = config->read<tuples3<double> >("D_bdry_vel");
    uD[0] = v.a; uD[1] = v.b; uD[2] = v.c;

    // external force (fluid)
    v = config->read<tuples3<double> >("ext_force_fluid");
    extF[0] = v.a; extF[1] = v.b; extF[2] = v.c;

    // heat transfer/thermal dynamics
    if(heat_xfer){
      th_diff_fl = config->read<double>("therm_diff_fl");
      Ra = config->read<double>("Rayleigh_number");
      tuples6<double> thetaBdry(-1.,-1.,-1.,-1.,-1.,-1.);
      thetaBdry = config->read<tuples6<double> >("bdry_temperature");
      bdryTemp[0] = thetaBdry.a; bdryTemp[1] = thetaBdry.b;
      bdryTemp[2] = thetaBdry.c; bdryTemp[3] = thetaBdry.d;
      bdryTemp[4] = thetaBdry.e; bdryTemp[5] = thetaBdry.f;
      T_hot = config->read<double>("T_hot");
      T_cold = config->read<double>("T_cold");
    }
  }
  // nu_b, temperature, fluct_seed
  if(lb_model=="mrt"){
    if(!plbpd->RSTRT){
      // nu_b
      nu_b = config->read<double>("nu_b");
      // whether or not fluctuations are present
      string answer = config->read<string>("fluctuations");
      if(answer=="yes")
	fluctuations = true;
      else if(answer=="no")
	fluctuations = false;
      else{
	string str = "input param 'fluctuations' can be one of yes/no";
	error->all(str);
      }
      if(fluctuations){
	// Boltzmann constant
	kB = config->read<double>("kB");
	// temperature
	temperature = config->read<double>("temperature");
	// seed for random fluctuations
	fluct_seed = config->read<int>("fluct_seed");
      }
    }
  }

  // nAg
  nAg = config->read<int>("num_aglmrts");
  if((nAg>0) and (!plbpd->RSTRT)){
    // is the aglmrt fixed?
    string answer = config->read<string>("aglmrt_fixed");
    if(answer=="yes")
      aglmrt_fixed = true;
    else if(answer=="no")
      aglmrt_fixed = false;
    else{
      string str = "input param 'aglmrt_fixed' can be one of yes/no";
      error->all(str);
    }
    if(!aglmrt_fixed){
      // integrator
      integrator = config->read<string>("integrator");
      // p_start_step
      p_start_step = config->read<int>("p_start_step");
      // external force on aglmrt
      v = config->read<tuples3<double> >("ext_force_aglmrt");
      extFaglmrt[0] = v.a; extFaglmrt[1] = v.b; extFaglmrt[2] = v.c;      
    }
    // aglmrt data file
    aggDataFile = config->read<string>("ag_data_file");
  }
}

/******************************************************************************/

void Input::read_restart_inp(){
  // read file name and open it
  rstrtDataFile = config->read<string>("rstrt_data_file");
  hid_t file_id = H5Fopen(rstrtDataFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id<0){
    string str = "file " + rstrtDataFile + " not found";
    error->all(str);
  }
  
  // read the parameters
  restart->read_inpparams(file_id);
  
  // close file
  H5Fclose(file_id);
}

/******************************************************************************/

void Input::check(){
  string DIRECTION[] = {"EAST", "WEST", "NORTH", "SOUTH", "UP", "DOWN"};

  if((Lx<=0) or (Ly<=0) or (Lz<=0)){
    string str("field 'domain_dim' not set correctly");
    cout<<Lx<<" "<<Ly<<" "<<Lz<<endl;
    error->all(str);
  }
  if((Lx==1) or (Lz==1)){
    string str("for 2D simulations, set Ly to 1, Lx and Lz greater than 1");
    error->all(str);
  }
  if((Ly==1) and bdryType[2]!="p"){
    string str("for 2D simulations, 'bdry_type' must be 'p' in the y-dir");
    error->all(str);
  }
  if((Ly==1) and (npy>1)){
    string str("for 2D simulations (Ly=1), npy can only be 1");
    error->all(str);
  }

  // check periodicity
  if( (bdryType[0]=="p") or
      (bdryType[1]=="p") ){
    if(bdryType[0]!=bdryType[1]){
      fin.close();
      string str("periodic and non-periodic bdries in the x direction");
      error->all(str);
    }
  } 
  if( (bdryType[2]=="p") or
      (bdryType[3]=="p") ){
    if(bdryType[2]!=bdryType[3]){
      fin.close();
      string str("periodic and non-periodic bdries in the y direction");
      error->all(str);
    }
  }
  if( (bdryType[4]=="p") or
      (bdryType[5]=="p") ){
    if(bdryType[4]!=bdryType[5]){
      fin.close();
      string str("periodic and non-periodic bdries in the z direction");
      error->all(str);
    }
  }
  // bdryType can be one of: i (inflow), o (outflow), n (noslip)
  if( (bdryType[0]!="o") and
      (bdryType[0]!="n") and
      (bdryType[1]!="p") ){
    fin.close();
    string str("east bdry can be one of o/n/p");
    error->all(str);
  }
  if( (bdryType[1]!="i") and
      (bdryType[1]!="n") and
      (bdryType[1]!="p") ){
    fin.close();
    string str("west bdry can be one of i/n/p");
    error->all(str);
  }
  for(int i=2; i<6; i++){
    if( (bdryType[i]!="n") and
	(bdryType[i]!="p") ){
      fin.close();
      string str(DIRECTION[i] + " bdry can be one of n/p");
      error->all(str);
    }
  }
  if( ((bdryType[0]=="o") and (bdryType[1]!="i")) or
      ((bdryType[1]=="i") and (bdryType[0]!="o")) ){
    fin.close();
    string str("for inflow/outflow: west/east bdry must be i/o");
    error->all(str);
  }

  // bdryType_th can be one of: d (dirichlet), n (neumann), p(periodic),
  //                            i (inflow), o (outflow)
  // for inflow/outlfow/periodic, no thermal BCs are applied
  // also, if fluid BCs are periodic in a particular direction, the thermal
  // BCs have to be periodic in that direction. so is the case for 
  // inflow/outflow
  // ONLY 'NOSLIP' BOUNDARY CAN BE DIRICHLET/NEUMANN
  if(heat_xfer){
    // TODO
    




  }

  if((lb_model!="bgk") and (lb_model!="mrt")){
    string str("input param 'lb_model' can be one of bgk/mrt");
    error->all(str);
  }

  if(rhoFluid<0.0){
    string str("field 'rho_fluid' not set correctly");
    error->all(str);
  }

  if(rhoSD<0.0){
    string str("field 'rho_solid_domain' not set correctly");
    error->all(str);
  }

  if(nu<0.0){
    string str("field 'nu' not set correctly");
    error->all(str);
  }

  if(th_diff_fl<0.0){
    string str("field 'th_diff_fl' not set correctly");
    error->all(str);
  }

  if(lb_model=="mrt"){
    // nu_b
    if(nu_b<0.0){
      string str("field 'nu_b' not set correctly");
      error->all(str);
    }
    // temperature
    if((fluctuations) and (temperature<=0.0)){
      string str("field 'temperature' not set correctly");
      error->all(str);
    }
  }

  if(nAg<0){
    string str("field 'nAg' not set correctly");
    error->all(str);
  }

  if((nAg>0) and (!aglmrt_fixed)){
    if((integrator!="rk4") and (integrator!="rkf45")){
      string str("input param 'integrator' can be one of rk4/rkf45");
      error->all(str);
    }

    if(p_start_step<1){
      string str("input param 'p_start_step' must be > 1");
      error->all(str);
    }
  }

  if(nSteps<=0){
    string str("field 'n_steps' not set correctly");
    error->all(str);
  }

  if(nWriteVel<=0){
    string str("field 'n_write_vel' not set correctly");
    error->all(str);
  }

  if(nWriteRestart<=0){
    string str("field 'n_write_state' not set correctly");
    error->all(str);
  }

  if(nPrintScreen<=0){
    string str("field 'n_print_screen' not set correctly");
    error->all(str);
  }

  if(nAg>0){
    if(nWriteForce<=0){
      string str("field 'n_write_force' not set correctly");
      error->all(str);
    }
    if(!aglmrt_fixed){
      if(nWriteXYZ<=0){
	string str("field 'n_write_xyz' not set correctly");
	error->all(str);
      }
      if(nWriteAVel<=0){
	string str("field 'n_write_a_vel' not set correctly");
	error->all(str);
      }
    }
  }

  if(heat_xfer){
    if(th_diff_fl<0){
      string str = "field 'th_diff_fl' cannot be negative";
      error->all(str);
    }
    if(T_hot<0){
      string str = "field 'T_hot' cannot be negative";
      error->all(str);
    }
    if(T_cold<0){
      string str = "field 'T_cold' cannot be negative";
      error->all(str);
    }
    if(T_hot<T_cold){
      string str = "oops! T_hot<T_cold";
      error->all(str);
    }
    if(Ra<0){
      string str = "field 'Rayleigh_number' cannot be negative";
      error->all(str);
    }
  }
}

/******************************************************************************/

void Input::read_restart_data(){
  // open the restart file
  hid_t file_id = H5Fopen(rstrtDataFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id<0){
    string str("file '" + rstrtDataFile + "' not found");
    error->all(str);
  }
  
  // read data
  restart->read_data(file_id);
  
  // close file
  H5Fclose(file_id);
}

/******************************************************************************/

void Input::read_steps(){
  // nSteps, nWriteVel, nWriteRestart, nPrintScreen, nWriteForce,
  // nWriteXYZ, nWriteAVel
  nSteps = config->read<int>("n_steps");
  nWriteVel = config->read<int>("n_write_vel");
  nWriteRestart = config->read<int>("n_write_rstrt");
  nPrintScreen = config->read<int>("n_print_screen");
  if(nAg>0){
    nWriteForce = config->read<int>("n_write_force");
    if(!aglmrt_fixed){
      nWriteXYZ = config->read<int>("n_write_xyz");
      nWriteAVel = config->read<int>("n_write_ag_vel");
    }
  }
}

/******************************************************************************/
