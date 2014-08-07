#include <iostream>
#include <cstring>

#include "mpi.h"
#include "plbpd.hpp"
#include "config.hpp"
#include "input.hpp"
#include "comm.hpp"
#include "error.hpp"
#include "memory.hpp"
#include "output.hpp"
#include "restart.hpp"
#include "timer.hpp"
#include "model.hpp"
#include "domain.hpp"
#include "boundary.hpp"
#include "lattice.hpp"
#include "lb_dynamics.hpp"
#include "bgk.hpp"
#include "mrt.hpp"
#include "th_dynamics.hpp"
#include "p_dynamics.hpp"

/******************************************************************************/

PLBPD::PLBPD(int argc, char **argv, MPI_Comm communicator, MPI_Info info){

  try{
    config = new Config(this);
    input = new Input(this);
    comm = new Comm(this, communicator, info); world = comm->world;
    error = new Error(this);
    memory = new Memory(this);
    output = new Output(this);
    restart = new Restart(this);
    timer = new Timer(this);
    model = new Model(this);
    domain = new Domain(this);
    boundary = new Boundary(this);
    lattice = new Lattice(this);
    p_dynamics = new P_Dynamics(this);
    th_dynamics = new TH_Dynamics(this);
    lb_dynamics = NULL;
  }
  catch(std::bad_alloc &xa){
    std::cerr<<"PLBPD Ctor failed to allocate memory\n";
  }

  // check input variables
  // check input options and get the input/restart file name
  if(argc<3){
    string str = "usage is\n      either ./plbpd -i <input file>";
    str += "\n      or     ./plbpd -r <restart file>";
    error->notify(str);
    str = "incorrect options";
    error->all(str);
  }
  if((strcmp(argv[1],"-i")!=0) and (strcmp(argv[1],"-r")!=0)){
    string str = "usage is\n      either ./plbpd -i <input file>";
    str += "\n      or     ./plbpd -r <restart file>";
    error->notify(str);
    str = "incorrect options";
    error->all(str);
  }
  infile.append(argv[2]);
  if(strcmp(argv[1],"-r")==0)
    RSTRT = true;
  else
    RSTRT = false;
}

/******************************************************************************/

PLBPD::~PLBPD(){
  if(lb_dynamics) delete lb_dynamics;
  if(th_dynamics) delete th_dynamics;
  if(p_dynamics) delete p_dynamics;
  if(lattice) delete lattice;
  if(boundary) delete boundary;
  if(domain) delete domain;
  if(model) delete model;
  if(timer) delete timer;
  if(restart) delete restart;
  if(output) delete output;
  if(memory) delete memory;
  if(error) delete error;
  if(comm) delete comm;
  if(input) delete input;
  if(config) delete config;
}

/******************************************************************************/

void PLBPD::setup(){
  /* DO NOT CHANGE THE ORDER OF CALLS TO SETUPS*/
  input->read();             // common read for both new and restart
  if (input->lb_model=="bgk")
    lb_dynamics = new BGK(this);
  else if (input->lb_model=="mrt")
    lb_dynamics = new MRT(this);
  else{
    string str = "LB model " + input->lb_model + " not recognized";
    error->all(str);
  }
  restart->setup();
  if(RSTRT)
    input->read_restart_inp();
  input->read_steps();
  input->check();

  timer->setup();
  model->setup();            // needs input->heat_xfer
  comm->setup();             // needs model->fluid/thermal->nVelocity
  domain->setup();           // input->Lx, needs comm->npx, comm->myloc[0]
  restart->alloc();          // needs domain->Lx
  comm->alloc();             // needs domain->lx etc.
  lattice->setup();          // needs domain->lx etc.
  boundary->setup();         // needs domain and lattice
  if(input->heat_xfer)
    boundary->setup_heatxfer();
  lb_dynamics->setup();
  if(input->heat_xfer)
    th_dynamics->setup();
  p_dynamics->setup();
  output->setup();
}
  
/******************************************************************************/

void PLBPD::run(){
  // some initializations
  if(RSTRT){
    input->read_restart_data();
  }
  else{
    timestep = 0;
    lb_dynamics->init();
    lb_dynamics->calc_moments();
    if(input->heat_xfer){
      th_dynamics->init();     // after lb_dynamics->init(), needs u @ node
      boundary->apply_thermalBC();
      th_dynamics->calc_temperature();
    }
    if(input->nAg>0){
      p_dynamics->init();
      p_dynamics->reset_flag_fluidist();
      p_dynamics->calc_mass_com_moi();
    }
  }

  comm->communicate();

  output->nodecount();
  output->header();
  output->print2screen(timestep);
  output->write_moments(timestep);
  if((input->nAg>0) and (!input->aglmrt_fixed)){
    output->write_prmry_xyz(timestep);
    output->write_aglmrt_vel(timestep);
  }

  output->write_restart(timestep);

  // main loop
  timer->start(MAINLOOP);
  for(tStep=timestep+1; tStep<timestep+input->nSteps+1; tStep++){
    /******LB-MD-heat transfer steps****/
    if(input->heat_xfer){
      th_dynamics->collideNstream();
      boundary->apply_thermalBC();
      th_dynamics->calc_temperature();
    }
    lb_dynamics->collideNstream();
    boundary->apply_fluidBC(); // hydrodynamic forces are calculated
                               // here when aglmrts are present

    if((input->nAg>0) and
       (!input->aglmrt_fixed) and
       (tStep>=input->p_start_step)){
      p_dynamics->integrate();
      p_dynamics->reset_flag_fluidist();
    }
    lb_dynamics->calc_moments();

    comm->communicate();
    /**end LB-MD-heat transfer steps****/

    /*******write output****************/
    timer->start(OUTPUT);
    if(tStep%input->nPrintScreen==0)
      output->print2screen(tStep);
    if(input->nAg>0){
      if(tStep%input->nWriteForce==0){
	output->write_force(tStep);
	output->write_torque(tStep);
      }
      if(!input->aglmrt_fixed){
	if(tStep%input->nWriteXYZ==0)
	  output->write_prmry_xyz(tStep);
	if(tStep%input->nWriteAVel==0)
	  output->write_aglmrt_vel(tStep);
      }
    }
    if(tStep%input->nWriteVel==0)
      output->write_moments(tStep);
    if(tStep%input->nWriteRestart==0)
      output->write_restart(tStep);

    timer->end(OUTPUT);
    /*****end write output**************/
  }
  timer->end(MAINLOOP);
 
  if((tStep-1)%input->nWriteVel!=0)
    output->write_moments(tStep-1);
  if((tStep-1)%input->nWriteRestart!=0)
    output->write_restart(tStep-1);  
  output->write_time_taken(); // writes elapsed times for each operation
}

/******************************************************************************/
