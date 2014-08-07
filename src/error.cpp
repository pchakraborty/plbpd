#include <cstdlib>
#include <iostream>

#include "mpi.h"
#include "error.hpp"
#include "comm.hpp"
#include "plbpd.hpp"

using std::cout;
using std::endl;

/******************************************************************************/

Error::Error(PLBPD *plbpd): Pointers(plbpd){
}

/******************************************************************************/

Error::~Error(){
}

/******************************************************************************/

void Error::all(string &errstr){
  if(comm->me==0) cout<<"ERROR: "<<errstr<<endl;
  delete plbpd;

  MPI_Finalize();
  exit(1);
}

/******************************************************************************/

void Error::warning(string &warnstr){
  if(comm->root) cout<<"WARNING: "<<warnstr<<endl;
}
  
/******************************************************************************/

void Error::notify(string &notestr){
  if(comm->root) cout<<"NOTE: "<<notestr<<endl;
}
  
/******************************************************************************/
