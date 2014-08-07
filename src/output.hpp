#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "pointers.hpp"
#include <fstream>
#include <hdf5.h>

class Output: protected Pointers{
public:
  Output(PLBPD* plbpd);
  ~Output();
  void setup();
  void header();
  void nodecount();
  void print2screen(int tStep);
  void write_moments(int tStep);
  void write_restart(int tStep);
  void write_force(int tStep);
  void write_torque(int tStep);
  void write_nodetype_xz();
  void write_nodetype_xy();
  void write_time_taken();
  void write_prmry_xyz(int tStep);
  void write_aglmrt_vel(int tStep);

  std::ofstream fout;    // output (rho,v)
  std::ofstream rout;    // output (state - binary output)
  std::ofstream Fout;    // output (hyd force)
  std::ofstream Tout;    // output (hyd torque)
  std::ofstream tout;    // output (aglmrt trajectory)
  std::ofstream vout;    // output (aglmrt velocity)
  std::ofstream boutxz;  // boundary
  
private:
  int nVelocity;
  std::string forceFile;
  std::string torqueFile;
  std::string prmryTrajFile;
  std::string aglmrtVelFile;
  FILE *logfile;

  /************ for writing binary data************************************/
  // data and memory spaces for type/rho, u, PI, n, input parameters
  hid_t dataspace,dataspace_3,dataspace_6,dataspace_nVel,dataspace_inp;
  hid_t memspace, memspace_3, memspace_6, memspace_nVel;
  hid_t inp_params_id;
};

#endif
