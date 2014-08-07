#ifndef RESTART_HPP
#define RESTART_HPP

#include "pointers.hpp"
#include "utils.hpp"
#include "global.hpp"
using namespace PLBPD_NS;

#include <hdf5.h>

class Restart: protected Pointers{
public:
  // structure to hold some input variables
  struct inp_params{
    int tstep;
    int npx, npy, npz;
    int Lx, Ly, Lz;
    char bdryType[6];
    char lb_model[10];
    double nu;
    double uE[3];
    double uW[3];
    double uN[3];
    double uS[3];
    double uU[3];
    double uD[3];
    double extF[3];
    // for mrt
    double nu_b;
    double temperature;
    int fluct_seed;
    // for aglmrts
    int nAg;
    int aglmrt_fixed;
    int p_start_step;
    char integrator[10];
    // for non-isothermal flow
    int heat_xfer;
    double Ra;
    double th_diff_fl;
    char bdryType_th[6];
    double bdryTemp[6];
  };

  struct ag_data{
    unsigned int ID;
    unsigned int nPr;
    // each primary consists of 7 doubles, r, c[3] and cb[3]
    double prmry[PLBPD_NS::MAXPRMRYS*PLBPD_NS::NPRMRYDBLES];
    double rho;
    double vol;
    double th_diff;
    double mass;
    double Ib[3];
    double IbI[3];
    double com[3];
    double R[9];
    double P[3];
    double L[3];
    double Iinv[9];
    double u[3];
    double omega[3];
    double F_hyd[3];
    double T_hyd[3];
    double F_tot[3];
    double T_tot[3];
  };

  Restart(PLBPD *);
  virtual ~Restart();

  void setup();
  void alloc();
  void read_inpparams(hid_t file_id);
  void read_data(hid_t file_id);
  void write_moments(hid_t file_id, hid_t plist_id);
  void write_restart(hid_t file_id, hid_t plist_id, int timestep);

private:
  // data and memory spaces for type/rho, u, n, input parameters
  hid_t filespace, filespace_3, filespace_nV_fl, filespace_nV_th;
  hid_t filespace_inp, filespace_ag;
  hid_t memspace, memspace_3, memspace_nV_fl, memspace_nV_th;
  hid_t inp_params_id;
  hid_t agdata_id;  
};

#endif
