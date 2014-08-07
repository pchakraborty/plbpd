#include "lb_dynamics.hpp"
#include "model.hpp"
#include "lattice.hpp"
#include "domain.hpp"
#include "input.hpp"
#include "comm.hpp"

#include "mpi.h"

/******************************************************************************/

LB_Dynamics::LB_Dynamics(PLBPD *plbpd): Pointers(plbpd){
};

/******************************************************************************/

LB_Dynamics::~LB_Dynamics(){
}

/******************************************************************************/
