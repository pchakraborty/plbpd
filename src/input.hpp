#ifndef INPUT_HPP
#define INPUT_HPP

#include "pointers.hpp"
#include <fstream>

class Input: protected Pointers{
public:
  Input(PLBPD *);
  ~Input();
  void read();
  void read_steps();
  void read_restart_inp();
  void read_restart_data();
  void check();
  
  // Comm parameters
  int npx, npy, npz;

  // Domain parameters
  int Lx, Ly, Lz;
  std::string bdryType[6];

  // LB_Dynamics parameters
  std::string lb_model;  // which model to use, BGK, MRT etc.
  double rhoFluid;       // fluid density
  double rhoSD;          // solid domain density
  double nu;             // viscosity
  double flowVel[3];     // flow velocity
  double uE[3];          // E bdry velocity
  double uW[3];          // W bdry velocity
  double uN[3];          // N bdry velocity
  double uS[3];          // S bdry velocity
  double uU[3];          // U bdry velocity
  double uD[3];          // D bdry velocity
  double extF[3];        // external force

  // for lb_model=="mrt"
  double nu_b;           // bulk viscosity
  double kB;             // boltzmann constant
  double temperature;    // fluid temperature (Kelvin)
  bool fluctuations;
  unsigned int fluct_seed; // seed for generating random numbers

  // heat transfer/thermal dynamics
  bool heat_xfer;        // whether or not to compute thermal diffusion
  double th_diff_fl;     // thermal diffusion coefficient of fluid
  double T_cold, T_hot;  // cold and hot temperatures
  std::string bdryType_th[6]; // dirichlet/neumann
  double bdryTemp[6];    // boundary temperature E, W, N, S, U, D
  double Ra;             // Rayleigh number
  
  // run and output parameters
  int nSteps;            // num steps of LB
  int nWriteVel;         // write rho, vel to file
  int nWriteRestart;     // write state of the system to file
  int nPrintScreen;      // write to screen
  int nWriteForce;       // write LB force on agglomerates to file
  int nWriteXYZ;         // write primary centre coords to xyz file
  int nWriteAVel;        // write aglmrt velocity data to file

  // P_Dynamics parameters
  std::string integrator;// can be rk4 or rkf45
  //int n_pdsteps;       // num of PD steps per LB step
  int nAg;               // num of agglomerates
  bool aglmrt_fixed;     // whether the aglmrt is fixed or moving
  int p_start_step;      // time step from when md starts
  double extFaglmrt[3];  // external force on aglmrt

  std::string aggDataFile;
  std::string rstrtDataFile;
  
  
  std::ifstream fin;
};

#endif
