#include "comm.hpp"
#include "error.hpp"
#include "domain.hpp"
#include "input.hpp"
#include "model.hpp"
#include "lattice.hpp"
#include "boundary.hpp"
#include "timer.hpp"
#include "utils.hpp"

#include <algorithm>

/******************************************************************************/

enum{EASTWEST, NORTHSOUTH, UPDOWN};

/******************************************************************************/

Comm::Comm(PLBPD *plbpd, MPI_Comm communicator, MPI_Info i): Pointers(plbpd){
  
  world = communicator;
  info = i;

  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);

  root = false;
  if(0==me) root = true;

  // Input::file() sets the values of npx, npy, npz
  npx = 0;
  npy = 0;
  npz = 0;

  sendbuf_f = NULL;
  sendbuf_b = NULL;
  recvbuf_f = NULL;
  recvbuf_b = NULL;
}

/******************************************************************************/

Comm::~Comm(){
  if(sendbuf_f) delete [] sendbuf_f;
  if(sendbuf_b) delete [] sendbuf_b;
  if(recvbuf_f) delete [] recvbuf_f;
  if(recvbuf_b) delete [] recvbuf_b;
}

/******************************************************************************/

void Comm::setup(){
  // complete MPI setup, after input file has been read
  npx = input->npx;
  npy = input->npy;
  npz = input->npz;

  int ncells = npx*npy*npz;

  if(nprocs!=ncells){
    string str = "num of procs "+utils::T2str<int>(nprocs)+" and num of cells "+
      utils::T2str<int>(ncells)+" do not match";
    error->all(str);
  }
  
  /* compute nbring proc numbers in the
     east/west, north/south and up/down directions */
  myloc[0] = me%npx;
  myloc[1] = (me/npx)%npy;
  myloc[2] = me/(npx*npy);

  // in the following, its input->bdryType[i] instead of
  // boundary->type[i], since boundary->type[i]
  // has not been set yet
  if(myloc[0]!=npx-1)
    east = me+1;
  else{
    if(input->bdryType[0]=="p"){
      east = 0 + npx * (myloc[1] + npy*myloc[2]);
    }
    else
      east = MPI_PROC_NULL; // no nbr in that direction
  }
  if(myloc[0]!=0)
    west = me-1;
  else{
    if(input->bdryType[1]=="p")
      west = (npx-1) + npx * (myloc[1] + npy*myloc[2]);
    else
      west = MPI_PROC_NULL;
  }
  if(myloc[1]!=npy-1)
    north = me+npx;
  else{
    if(input->bdryType[2]=="p")
      north = myloc[0] + npx * (0 + npy*myloc[2]);
    else
      north = MPI_PROC_NULL;
  }
  if(myloc[1]!= 0)
    south = me-npx;
  else{
    if(input->bdryType[3]=="p")
      south = myloc[0] + npx *((npy-1) + npy*myloc[2]);
    else
      south = MPI_PROC_NULL;
  }
  if(myloc[2]!=npz-1)
    up = me+npx*npy;
  else{
    if(input->bdryType[4]=="p")
      up = myloc[0] + npx * (myloc[1] + npy*0);
    else
      up = MPI_PROC_NULL;
  }
  if(myloc[2]!= 0)
    down = me-npx*npy;
  else{
    if(input->bdryType[5]=="p")
      down = myloc[0] + npx * (myloc[1] + npy*(npz-1));
    else
      down = MPI_PROC_NULL;
  }

  /* 
     set nDoubles: num of 'doubles' at each node that need to be
     communicated,
     type, rho, u[3], n[model->nVelocity]

     actually, we do not need to pass all the nk's. only the unknowns
     at the receiving bdry. so it can be further optimized

     during the first call to communicate, however, all the nk's
     need to be communicated

     for heat transfer, we need to communicate some more data
     u[3] + 2*(temperature, alpha, nth[model->thermal->nVelocity])
  */
  nDoubles = 1 + 1 + 3 + model->fluid->nVelocity;
  if(input->heat_xfer)
    nDoubles += 3 + 2*(1 + 1 + model->thermal->nVelocity);
  // additional memory for random variables in case of mrt
  if((input->lb_model=="mrt") and (input->fluctuations))
    nDoubles += model->fluid->nVelocity;
}

/******************************************************************************/

void Comm::alloc(){
  int maxlx, maxly, maxlz;
  MPI_Allreduce(&domain->lx,&maxlx,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&domain->ly,&maxly,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&domain->lz,&maxlz,1,MPI_INT,MPI_MAX,world);

  int size = std::max(maxlx,maxly);
  size = std::max(maxlz,size);

  // +4 accounts for the boundary layers
  // +3 accounts for xdim, ydim, zdim
  maxbufsize = (size+4)*(size+4)*nDoubles + 3;

  try{
    sendbuf_f = new double[maxbufsize];
    sendbuf_b = new double[maxbufsize];
    recvbuf_f = new double[maxbufsize];
    recvbuf_b = new double[maxbufsize];
  }
  catch(std::bad_alloc &xa){
    string str("Comm::alloc() failed to allocate memory");
    error->all(str);
  }
}

/******************************************************************************/

void Comm::pack(const int direction){
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  int nV_fl = model->fluid->nVelocity;
  int nV_th = 0;
  if(input->heat_xfer) 
    nV_th = model->thermal->nVelocity;

  int ***nodetype = lattice->type;
  double ****n = lattice->n;
  double ***rho = lattice->rho;
  double ****u = lattice->u;
  double ***theta = lattice->theta;
  double ****nth = lattice->nth;
  double ***alpha = lattice->alpha;
  double ****RV = lattice->RV;

  int xl, yl, zl, ctr;

  if(direction==EASTWEST){
    // east buffer
    sendbuf_f[0] = static_cast<double>(xdim);
    sendbuf_f[1] = static_cast<double>(ydim);
    sendbuf_f[2] = static_cast<double>(zdim);
    ctr = 3;

    xl = xdim + 1;
    // for(zl=1; zl<zdim+3; zl++){
    //   for(yl=1; yl<ydim+3; yl++){
    for(zl=0; zl<zdim+4; zl++){
      for(yl=0; yl<ydim+4; yl++){
	sendbuf_f[ctr++] = static_cast<double>(nodetype[zl][yl][xl]);
	sendbuf_f[ctr++] = rho[zl][yl][xl];
	for(int i=0; i<3; i++)
	  sendbuf_f[ctr++] = u[zl][yl][xl][i];
	for(int k=0; k<nV_fl; k++)
	  sendbuf_f[ctr++] = n[zl][yl][xl][k];
	// random numbers for mrt
	if((input->lb_model=="mrt") and (input->fluctuations)){
	  for(int i=0; i<nV_fl; i++)
	    sendbuf_f[ctr++] = RV[zl][yl][xl][i];
	}
	// heat transfer
	if(input->heat_xfer){
	  for(int i=0; i<3; i++)
	    sendbuf_f[ctr++] = u[zl][yl][xl-1][i];
	  sendbuf_f[ctr++] = theta[zl][yl][xl];
	  sendbuf_f[ctr++] = theta[zl][yl][xl-1];
	  sendbuf_f[ctr++] = alpha[zl][yl][xl];
	  sendbuf_f[ctr++] = alpha[zl][yl][xl-1];
	  for(int k=0; k<nV_th; k++){
	    sendbuf_f[ctr++] = nth[zl+1][yl+1][xl+1][k];
	    sendbuf_f[ctr++] = nth[zl+1][yl+1][xl+1-1][k];
	  }
	}
      }
    }
    count_f = ctr;
  
    // west buffer
    sendbuf_b[0] = static_cast<double>(xdim);
    sendbuf_b[1] = static_cast<double>(ydim);
    sendbuf_b[2] = static_cast<double>(zdim);
    ctr = 3;  

    xl = 2;
    // for(zl=1; zl<zdim+3; zl++){
    //   for(yl=1; yl<ydim+3; yl++){
    for(zl=0; zl<zdim+4; zl++){
      for(yl=0; yl<ydim+4; yl++){
	sendbuf_b[ctr++] = static_cast<double>(nodetype[zl][yl][xl]);
	sendbuf_b[ctr++] = rho[zl][yl][xl];
	for(int i=0; i<3; i++)
	  sendbuf_b[ctr++] = u[zl][yl][xl][i];
	for(int k=0; k<nV_fl; k++)
	  sendbuf_b[ctr++] = n[zl][yl][xl][k];
	// random numbers for mrt
	if((input->lb_model=="mrt") and (input->fluctuations)){
	  for(int i=0; i<nV_fl; i++)
	    sendbuf_b[ctr++] = RV[zl][yl][xl][i];
	}
	// heat transfer
	if(input->heat_xfer){
	  for(int i=0; i<3; i++)
	    sendbuf_b[ctr++] = u[zl][yl][xl+1][i];
	  sendbuf_b[ctr++] = theta[zl][yl][xl];
	  sendbuf_b[ctr++] = theta[zl][yl][xl+1];
	  sendbuf_b[ctr++] = alpha[zl][yl][xl];
	  sendbuf_b[ctr++] = alpha[zl][yl][xl+1];
	  for(int k=0; k<nV_th; k++){
	    sendbuf_b[ctr++] = nth[zl+1][yl+1][xl+1][k];
	    sendbuf_b[ctr++] = nth[zl+1][yl+1][xl+1+1][k];
	  }
	}
      }
    }
    count_b = ctr;
  }

  else if(direction==NORTHSOUTH){
    // north buffer
    sendbuf_f[0] = static_cast<double>(xdim);
    sendbuf_f[1] = static_cast<double>(ydim);
    sendbuf_f[2] = static_cast<double>(zdim);
    ctr = 3;

    yl = ydim + 1;
    // for(zl=1; zl<zdim+3; zl++){
    //   for(xl=1; xl<xdim+3; xl++){
    for(zl=0; zl<zdim+4; zl++){
      for(xl=0; xl<xdim+4; xl++){
	sendbuf_f[ctr++] = static_cast<double>(nodetype[zl][yl][xl]);
	sendbuf_f[ctr++] = rho[zl][yl][xl];
	for(int i=0; i<3; i++)
	  sendbuf_f[ctr++] = u[zl][yl][xl][i];
	for(int k=0; k<nV_fl; k++)
	  sendbuf_f[ctr++] = n[zl][yl][xl][k];
	// random numbers for mrt
	if((input->lb_model=="mrt") and (input->fluctuations)){
	  for(int i=0; i<nV_fl; i++)
	    sendbuf_f[ctr++] = RV[zl][yl][xl][i];
	}
	// heat transfer
	if(input->heat_xfer){
	  for(int i=0; i<3; i++)
	    sendbuf_f[ctr++] = u[zl][yl-1][xl][i];
	  sendbuf_f[ctr++] = theta[zl][yl][xl];
	  sendbuf_f[ctr++] = theta[zl][yl-1][xl];
	  sendbuf_f[ctr++] = alpha[zl][yl][xl];
	  sendbuf_f[ctr++] = alpha[zl][yl-1][xl];
	  for(int k=0; k<nV_th; k++){
	    sendbuf_f[ctr++] = nth[zl+1][yl+1][xl+1][k];
	    sendbuf_f[ctr++] = nth[zl+1][yl+1-1][xl+1][k];
	  }
	}
      }
    }
    count_f = ctr;
  
    // south buffer
    sendbuf_b[0] = static_cast<double>(xdim);
    sendbuf_b[1] = static_cast<double>(ydim);
    sendbuf_b[2] = static_cast<double>(zdim);
    ctr = 3;  

    yl = 2;
    // for(zl=1; zl<zdim+3; zl++){
    //   for(xl=1; xl<xdim+3; xl++){
    for(zl=0; zl<zdim+4; zl++){
      for(xl=0; xl<xdim+4; xl++){
	sendbuf_b[ctr++] = static_cast<double>(nodetype[zl][yl][xl]);
	sendbuf_b[ctr++] = rho[zl][yl][xl];
	for(int i=0; i<3; i++)
	  sendbuf_b[ctr++] = u[zl][yl][xl][i];
	for(int k=0; k<nV_fl; k++)
	  sendbuf_b[ctr++] = n[zl][yl][xl][k];
	// random numbers for mrt
	if((input->lb_model=="mrt") and (input->fluctuations)){
	  for(int i=0; i<nV_fl; i++)
	    sendbuf_b[ctr++] = RV[zl][yl][xl][i];
	}
	// heat transfer
	if(input->heat_xfer){
	  for(int i=0; i<3; i++)
	    sendbuf_b[ctr++] = u[zl][yl+1][xl][i];
	  sendbuf_b[ctr++] = theta[zl][yl][xl];
	  sendbuf_b[ctr++] = theta[zl][yl+1][xl];
	  sendbuf_b[ctr++] = alpha[zl][yl][xl];
	  sendbuf_b[ctr++] = alpha[zl][yl+1][xl];
	  for(int k=0; k<nV_th; k++){
	    sendbuf_b[ctr++] = nth[zl+1][yl+1][xl+1][k];
	    sendbuf_b[ctr++] = nth[zl+1][yl+1+1][xl+1][k];
	  }
	}
      }
    }
    count_b = ctr;
  }

  else if(direction==UPDOWN){
    // up buffer
    sendbuf_f[0] = static_cast<double>(xdim);
    sendbuf_f[1] = static_cast<double>(ydim);
    sendbuf_f[2] = static_cast<double>(zdim);
    ctr = 3;

    zl = zdim + 1;
    // for(yl=1; yl<ydim+3; yl++){
    //   for(xl=1; xl<xdim+3; xl++){
    for(yl=0; yl<ydim+4; yl++){
      for(xl=0; xl<xdim+4; xl++){
	sendbuf_f[ctr++] = static_cast<double>(nodetype[zl][yl][xl]);
	sendbuf_f[ctr++] = rho[zl][yl][xl];
	for(int i=0; i<3; i++)
	  sendbuf_f[ctr++] = u[zl][yl][xl][i];
	for(int k=0; k<nV_fl; k++)
	  sendbuf_f[ctr++] = n[zl][yl][xl][k];
	// random numbers for mrt
	if((input->lb_model=="mrt") and (input->fluctuations)){
	  for(int i=0; i<nV_fl; i++)
	    sendbuf_f[ctr++] = RV[zl][yl][xl][i];
	}
	// heat transfer
	if(input->heat_xfer){
	  for(int i=0; i<3; i++)
	    sendbuf_f[ctr++] = u[zl-1][yl][xl][i];
	  sendbuf_f[ctr++] = theta[zl][yl][xl];
	  sendbuf_f[ctr++] = theta[zl-1][yl][xl];
	  sendbuf_f[ctr++] = alpha[zl][yl][xl];
	  sendbuf_f[ctr++] = alpha[zl-1][yl][xl];
	  for(int k=0; k<nV_th; k++){
	    sendbuf_f[ctr++] = nth[zl+1][yl+1][xl+1][k];
	    sendbuf_f[ctr++] = nth[zl+1-1][yl+1][xl+1][k];
	  }
	}
      }
    }
    count_f = ctr;
  
    // down buffer
    sendbuf_b[0] = static_cast<double>(xdim);
    sendbuf_b[1] = static_cast<double>(ydim);
    sendbuf_b[2] = static_cast<double>(zdim);
    ctr = 3;  

    zl = 2;
    // for(yl=1; yl<ydim+3; yl++){
    //   for(xl=1; xl<xdim+3; xl++){
    for(yl=0; yl<ydim+4; yl++){
      for(xl=0; xl<xdim+4; xl++){
	sendbuf_b[ctr++] = static_cast<double>(nodetype[zl][yl][xl]);
	sendbuf_b[ctr++] = rho[zl][yl][xl];
	for(int i=0; i<3; i++)
	  sendbuf_b[ctr++] = u[zl][yl][xl][i];
	for(int k=0; k<nV_fl; k++)
	  sendbuf_b[ctr++] = n[zl][yl][xl][k];
	// random numbers for mrt
	if((input->lb_model=="mrt") and (input->fluctuations)){
	  for(int i=0; i<nV_fl; i++)
	    sendbuf_b[ctr++] = RV[zl][yl][xl][i];
	}
	// heat transfer
	if(input->heat_xfer){
	  for(int i=0; i<3; i++)
	    sendbuf_b[ctr++] = u[zl+1][yl][xl][i];
	  sendbuf_b[ctr++] = theta[zl][yl][xl];
	  sendbuf_b[ctr++] = theta[zl+1][yl][xl];
	  sendbuf_b[ctr++] = alpha[zl][yl][xl];
	  sendbuf_b[ctr++] = alpha[zl+1][yl][xl];
	  for(int k=0; k<nV_th; k++){
	    sendbuf_b[ctr++] = nth[zl+1][yl+1][xl+1][k];
	    sendbuf_b[ctr++] = nth[zl+1+1][yl+1][xl+1][k];
	  }
	}
      }
    }
    count_b = ctr;
  }

  else{
    string str="invalid communication direction "+utils::T2str<int>(direction);
    error->all(str);
  }
}

/******************************************************************************/

void Comm::unpack(const int direction){

  int xdim, ydim, zdim;
  int nV_fl = model->fluid->nVelocity;
  int nV_th = 0;
  if(input->heat_xfer)
    nV_th = model->thermal->nVelocity;

  int ***nodetype = lattice->type;
  double ****n = lattice->n;
  double ***rho = lattice->rho;
  double ****u = lattice->u;
  double ***theta = lattice->theta;
  double ****nth = lattice->nth;
  double ***alpha = lattice->alpha;
  double ****RV = lattice->RV;

  int xl, yl, zl;
  int ctr;

  if(direction==EASTWEST){
    // received from west
    if(west!=MPI_PROC_NULL){
      xdim = static_cast<int>(recvbuf_f[0]);
      ydim = static_cast<int>(recvbuf_f[1]);
      zdim = static_cast<int>(recvbuf_f[2]);
      ctr = 3;

      xl = 1;	
      // for(zl=1; zl<zdim+3; zl++){
      // 	for(yl=1; yl<ydim+3; yl++){
      for(zl=0; zl<zdim+4; zl++){
	for(yl=0; yl<ydim+4; yl++){
	  nodetype[zl][yl][xl] = static_cast<int>(recvbuf_f[ctr++]);
	  rho[zl][yl][xl] = recvbuf_f[ctr++];
	  for(int i=0; i<3; i++)
	    u[zl][yl][xl][i] = recvbuf_f[ctr++];
	  for(int k=0; k<nV_fl; k++)
	    n[zl][yl][xl][k] = recvbuf_f[ctr++];
	  // random numbers for mrt
	  if((input->lb_model=="mrt") and (input->fluctuations)){
	    for(int i=0; i<nV_fl; i++)
	      RV[zl][yl][xl][i] = recvbuf_f[ctr++];
	  }
	  // heat transfer
	  if(input->heat_xfer){
	    for(int i=0; i<3; i++)
	      u[zl][yl][xl-1][i] = recvbuf_f[ctr++];
	    theta[zl][yl][xl] = recvbuf_f[ctr++];
	    theta[zl][yl][xl-1] = recvbuf_f[ctr++];
	    alpha[zl][yl][xl] = recvbuf_f[ctr++];
	    alpha[zl][yl][xl-1] = recvbuf_f[ctr++];
	    for(int k=0; k<nV_th; k++){
	      nth[zl+1][yl+1][xl+1][k] = recvbuf_f[ctr++];
	      nth[zl+1][yl+1][xl+1-1][k] = recvbuf_f[ctr++];
	    }
	  }
	}
      }
    }

    // received from east
    if(east!=MPI_PROC_NULL){
      xdim = static_cast<int>(recvbuf_b[0]);
      ydim = static_cast<int>(recvbuf_b[1]);
      zdim = static_cast<int>(recvbuf_b[2]);
      ctr = 3;

      xl = domain->lx+2;
      // for(zl=1; zl<zdim+3; zl++){
      // 	for(yl=1; yl<ydim+3; yl++){
      for(zl=0; zl<zdim+4; zl++){
	for(yl=0; yl<ydim+4; yl++){
	  nodetype[zl][yl][xl] = static_cast<int>(recvbuf_b[ctr++]);
	  rho[zl][yl][xl] = recvbuf_b[ctr++];
	  for(int i=0; i<3; i++)
	    u[zl][yl][xl][i] = recvbuf_b[ctr++];
	  for(int k=0; k<nV_fl; k++)
	    n[zl][yl][xl][k] = recvbuf_b[ctr++];
	  // random numbers for mrt
	  if((input->lb_model=="mrt") and (input->fluctuations)){
	    for(int i=0; i<nV_fl; i++)
	      RV[zl][yl][xl][i] = recvbuf_b[ctr++];
	  }
	  // heat transfer
	  if(input->heat_xfer){
	    for(int i=0; i<3; i++)
	      u[zl][yl][xl+1][i] = recvbuf_b[ctr++];
	    theta[zl][yl][xl] = recvbuf_b[ctr++];
	    theta[zl][yl][xl+1] = recvbuf_b[ctr++];
	    alpha[zl][yl][xl] = recvbuf_b[ctr++];
	    alpha[zl][yl][xl+1] = recvbuf_b[ctr++];
	    for(int k=0; k<nV_th; k++){
	      nth[zl+1][yl+1][xl+1][k] = recvbuf_b[ctr++];
	      nth[zl+1][yl+1][xl+1+1][k] = recvbuf_b[ctr++];
	    }
	  }
	}
      }
    }
  }

  else if(direction==NORTHSOUTH){
    // received from south
    if(south!=MPI_PROC_NULL){
      xdim = static_cast<int>(recvbuf_f[0]);
      ydim = static_cast<int>(recvbuf_f[1]);
      zdim = static_cast<int>(recvbuf_f[2]);
      ctr = 3;

      yl = 1;
      // for(zl=1; zl<zdim+3; zl++){
      // 	for(xl=1; xl<xdim+3; xl++){
      for(zl=0; zl<zdim+4; zl++){
	for(xl=0; xl<xdim+4; xl++){
	  nodetype[zl][yl][xl] = static_cast<int>(recvbuf_f[ctr++]);
	  rho[zl][yl][xl] = recvbuf_f[ctr++];
	  for(int i=0; i<3; i++)
	    u[zl][yl][xl][i] = recvbuf_f[ctr++];
	  for(int k=0; k<nV_fl; k++)
	    n[zl][yl][xl][k] = recvbuf_f[ctr++];
	  // random numbers for mrt
	  if((input->lb_model=="mrt") and (input->fluctuations)){
	    for(int i=0; i<nV_fl; i++)
	      RV[zl][yl][xl][i] = recvbuf_f[ctr++];
	  }
	  // heat transfer
	  if(input->heat_xfer){
	    for(int i=0; i<3; i++)
	      u[zl][yl-1][xl][i] = recvbuf_f[ctr++];
	    theta[zl][yl][xl] = recvbuf_f[ctr++];
	    theta[zl][yl-1][xl] = recvbuf_f[ctr++];
	    alpha[zl][yl][xl] = recvbuf_f[ctr++];
	    alpha[zl][yl-1][xl] = recvbuf_f[ctr++];
	    for(int k=0; k<nV_th; k++){
	      nth[zl+1][yl+1][xl+1][k] = recvbuf_f[ctr++];
	      nth[zl+1][yl+1-1][xl+1][k] = recvbuf_f[ctr++];
	    }
	  }
	}
      }
    }
  
    // received from north
    if(north!=MPI_PROC_NULL){
      xdim = static_cast<int>(recvbuf_b[0]);
      ydim = static_cast<int>(recvbuf_b[1]);
      zdim = static_cast<int>(recvbuf_b[2]);
      ctr = 3;

      yl = domain->ly+2;
      // for(zl=1; zl<zdim+3; zl++){
      // 	for(xl=1; xl<xdim+3; xl++){
      for(zl=0; zl<zdim+4; zl++){
	for(xl=0; xl<xdim+4; xl++){
	  nodetype[zl][yl][xl] = static_cast<int>(recvbuf_b[ctr++]);
	  rho[zl][yl][xl] = recvbuf_b[ctr++];
	  for(int i=0; i<3; i++)
	    u[zl][yl][xl][i] = recvbuf_b[ctr++];
	  for(int k=0; k<nV_fl; k++)
	    n[zl][yl][xl][k] = recvbuf_b[ctr++];
	  // random numbers for mrt
	  if((input->lb_model=="mrt") and (input->fluctuations)){
	    for(int i=0; i<nV_fl; i++)
	      RV[zl][yl][xl][i] = recvbuf_b[ctr++];
	  }
	  // heat transfer
	  if(input->heat_xfer){
	    for(int i=0; i<3; i++)
	      u[zl][yl+1][xl][i] = recvbuf_b[ctr++];
	    theta[zl][yl][xl] = recvbuf_b[ctr++];
	    theta[zl][yl+1][xl] = recvbuf_b[ctr++];
	    alpha[zl][yl][xl] = recvbuf_b[ctr++];
	    alpha[zl][yl+1][xl] = recvbuf_b[ctr++];
	    for(int k=0; k<nV_th; k++){
	      nth[zl+1][yl+1][xl+1][k] = recvbuf_b[ctr++];
	      nth[zl+1][yl+1+1][xl+1][k] = recvbuf_b[ctr++];
	    }
	  }
	}
      }
    }
  }

  else if(direction==UPDOWN){
    // received from down
    if(down!=MPI_PROC_NULL){
      xdim = static_cast<int>(recvbuf_f[0]);
      ydim = static_cast<int>(recvbuf_f[1]);
      zdim = static_cast<int>(recvbuf_f[2]);
      ctr = 3;

      zl = 1;
      // for(yl=1; yl<ydim+3; yl++){
      // 	for(xl=1; xl<xdim+3; xl++){
      for(yl=0; yl<ydim+4; yl++){
	for(xl=0; xl<xdim+4; xl++){
	  nodetype[zl][yl][xl] = static_cast<int>(recvbuf_f[ctr++]);
	  rho[zl][yl][xl] = recvbuf_f[ctr++];
	  for(int i=0; i<3; i++)
	    u[zl][yl][xl][i] = recvbuf_f[ctr++];
	  for(int k=0; k<nV_fl; k++)
	    n[zl][yl][xl][k] = recvbuf_f[ctr++];
	  // random numbers for mrt
	  if((input->lb_model=="mrt") and (input->fluctuations)){
	    for(int i=0; i<nV_fl; i++)
	      RV[zl][yl][xl][i] = recvbuf_f[ctr++];
	  }
	  // heat transfer
	  if(input->heat_xfer){
	    for(int i=0; i<3; i++)
	      u[zl-1][yl][xl][i] = recvbuf_f[ctr++];
	    theta[zl][yl][xl] = recvbuf_f[ctr++];
	    theta[zl-1][yl][xl] = recvbuf_f[ctr++];
	    alpha[zl][yl][xl] = recvbuf_f[ctr++];
	    alpha[zl-1][yl][xl] = recvbuf_f[ctr++];
	    for(int k=0; k<nV_th; k++){
	      nth[zl+1][yl+1][xl+1][k] = recvbuf_f[ctr++];
	      nth[zl+1-1][yl+1][xl+1][k] = recvbuf_f[ctr++];
	    }
	  }
	}
      }
    }

    // received from up
    if(up!=MPI_PROC_NULL){
      xdim = static_cast<int>(recvbuf_b[0]);
      ydim = static_cast<int>(recvbuf_b[1]);
      zdim = static_cast<int>(recvbuf_b[2]);
      ctr = 3;

      zl = domain->lz+2;
      // for(yl=1; yl<ydim+3; yl++){
      // 	for(xl=1; xl<xdim+3; xl++){
      for(yl=0; yl<ydim+4; yl++){
	for(xl=0; xl<xdim+4; xl++){
	  nodetype[zl][yl][xl] = static_cast<int>(recvbuf_b[ctr++]);
	  rho[zl][yl][xl] = recvbuf_b[ctr++];
	  for(int i=0; i<3; i++)
	    u[zl][yl][xl][i] = recvbuf_b[ctr++];
	  for(int k=0; k<nV_fl; k++)
	    n[zl][yl][xl][k] = recvbuf_b[ctr++];
	  // random numbers for mrt
	  if((input->lb_model=="mrt") and (input->fluctuations)){
	    for(int i=0; i<nV_fl; i++)
	      RV[zl][yl][xl][i] = recvbuf_b[ctr++];
	  }
	  // heat transfer
	  if(input->heat_xfer){
	    for(int i=0; i<3; i++)
	      u[zl+1][yl][xl][i] = recvbuf_b[ctr++];
	    theta[zl][yl][xl] = recvbuf_b[ctr++];
	    theta[zl+1][yl][xl] = recvbuf_b[ctr++];
	    alpha[zl][yl][xl] = recvbuf_b[ctr++];
	    alpha[zl+1][yl][xl] = recvbuf_b[ctr++];
	    for(int k=0; k<nV_th; k++){
	      nth[zl+1][yl+1][xl+1][k] = recvbuf_b[ctr++];
	      nth[zl+1+1][yl+1][xl+1][k] = recvbuf_b[ctr++];
	    }
	  }
	}
      }
    }
  }
  
  else{
    string str="invalid communication direction "+utils::T2str<int>(direction);
    error->all(str);
  }
}

/******************************************************************************/
  
void Comm::communicate(){
  timer->start(COMM);

  MPI_Status status;
  int tag = 999;

  // EW communication
  {
    pack(EASTWEST);
    if(npx==1){
      if(boundary->type[0]=="p"){
	// copy send buffers to their corresponding recv buffers
	for(int i=0; i<count_f;i++)
	  recvbuf_f[i] = sendbuf_f[i];
	for(int i=0; i<count_b; i++)
	  recvbuf_b[i] = sendbuf_b[i];
      }
      else{
	//printf("EW-do nothing\n");
      }
    }
    else{
      if(myloc[0]%2==0){
	MPI_Send(sendbuf_f,count_f,   MPI_DOUBLE,east,tag,world);
	MPI_Recv(recvbuf_f,maxbufsize,MPI_DOUBLE,west,tag,world,&status);
	MPI_Send(sendbuf_b,count_b,   MPI_DOUBLE,west,tag,world);
	MPI_Recv(recvbuf_b,maxbufsize,MPI_DOUBLE,east,tag,world,&status);
      }
      else{
	MPI_Recv(recvbuf_f,maxbufsize,MPI_DOUBLE,west,tag,world,&status);
	MPI_Send(sendbuf_f,count_f,   MPI_DOUBLE,east,tag,world);
	MPI_Recv(recvbuf_b,maxbufsize,MPI_DOUBLE,east,tag,world,&status);
	MPI_Send(sendbuf_b,count_b,   MPI_DOUBLE,west,tag,world);
      }
    }
    unpack(EASTWEST);
    //printf("proc: %d, EW\n", me);
  }

  // NS communication
  {
    pack(NORTHSOUTH);
    if(npy==1){
      if(boundary->type[2]=="p"){
	// copy send buffers to their corresponding recv buffers
	for(int i=0; i<count_f; i++)
	  recvbuf_f[i] = sendbuf_f[i];
	for(int i=0; i<count_b; i++)
	  recvbuf_b[i] = sendbuf_b[i];
      }
      else{
	//printf("NS-do nothing\n");
      }
    }
    else{
      if(myloc[1]%2==0){
	MPI_Send(sendbuf_f,count_f,   MPI_DOUBLE,north,tag,world);
	MPI_Recv(recvbuf_f,maxbufsize,MPI_DOUBLE,south,tag,world,&status);
	MPI_Send(sendbuf_b,count_b,   MPI_DOUBLE,south,tag,world);
	MPI_Recv(recvbuf_b,maxbufsize,MPI_DOUBLE,north,tag,world,&status);
      }
      else{
	MPI_Recv(recvbuf_f,maxbufsize,MPI_DOUBLE,south,tag,world,&status);
	MPI_Send(sendbuf_f,count_f,   MPI_DOUBLE,north,tag,world);
	MPI_Recv(recvbuf_b,maxbufsize,MPI_DOUBLE,north,tag,world,&status);
	MPI_Send(sendbuf_b,count_b,   MPI_DOUBLE,south,tag,world);
      }
    }
    unpack(NORTHSOUTH);
  }

  // UD communication
  {
    pack(UPDOWN);
    if(npz==1){
      if(boundary->type[4]=="p"){
	// copy send buffers to their corresponding recv buffers
	for(int i=0; i<count_f; i++)
	  recvbuf_f[i] = sendbuf_f[i];
	for(int i=0; i<count_b; i++)
	  recvbuf_b[i] = sendbuf_b[i];
      }
      else{
	//printf("UD-do nothing\n");
      }
    }
    else{
      if(myloc[2]%2==0){
	MPI_Send(sendbuf_f,count_f,   MPI_DOUBLE,up,  tag,world);
	MPI_Recv(recvbuf_f,maxbufsize,MPI_DOUBLE,down,tag,world,&status);
	MPI_Send(sendbuf_b,count_b,   MPI_DOUBLE,down,tag,world);
	MPI_Recv(recvbuf_b,maxbufsize,MPI_DOUBLE,up,  tag,world,&status);
      }
      else{
	MPI_Recv(recvbuf_f,maxbufsize,MPI_DOUBLE,down,tag,world,&status);
	MPI_Send(sendbuf_f,count_f,   MPI_DOUBLE,up,  tag,world);
	MPI_Recv(recvbuf_b,maxbufsize,MPI_DOUBLE,up,  tag,world,&status);
	MPI_Send(sendbuf_b,count_b,   MPI_DOUBLE,down,tag,world);
      }
    }
    unpack(UPDOWN);
  }
  timer->end(COMM);
}

/******************************************************************************/

void Comm::communicate_serial(){
  int nV_fl = model->fluid->nVelocity;
  if(nprocs!=1){
    string str("communicate_serial() works only in the serial case");
    error->all(str);
  }

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;
  double ****n = lattice->n;
  int ***type = lattice->type;

  // EW
  for(int zl=1; zl<zdim+3; zl++){
    for(int yl=1; yl<ydim+3; yl++){
      type[zl][yl][1] = type[zl][yl][xdim+1];
      type[zl][yl][xdim+2] = type[zl][yl][2];
      for(int k=0; k<nV_fl; k++){
	n[zl][yl][1][k] = n[zl][yl][xdim+1][k];
	n[zl][yl][xdim+2][k] = n[zl][yl][2][k];
      }
    }
  }

  // NS
  for(int zl=1; zl<zdim+3; zl++){
    for(int xl=1; xl<xdim+3; xl++){
      type[zl][1][xl] = type[zl][ydim+1][xl];
      type[zl][ydim+2][xl] = type[zl][2][xl];
      for(int k=0; k<nV_fl; k++){
	n[zl][1][xl][k] = n[zl][ydim+1][xl][k];
	n[zl][ydim+2][xl][k] = n[zl][2][xl][k];
      }
    }
  }
  
  // UD
  for(int yl=1; yl<ydim+3; yl++){
    for(int xl=1; xl<xdim+3; xl++){
      type[1][yl][xl] = type[zdim+1][yl][xl];
      type[zdim+2][yl][xl] = type[2][yl][xl];
      for(int k=0; k<nV_fl; k++){
	n[1][yl][xl][k] = n[zdim+1][yl][xl][k];
	n[zdim+2][yl][xl][k] = n[2][yl][xl][k];
      }
    }
  }
}

/******************************************************************************/
