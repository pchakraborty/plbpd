#include "boundary.hpp"
#include "input.hpp"
#include "lattice.hpp"
#include "model.hpp"
#include "domain.hpp"
#include "error.hpp"
#include "lb_dynamics.hpp"
#include "th_dynamics.hpp"
#include "p_dynamics.hpp"
#include "timer.hpp"
#include "comm.hpp"
#include "utils.hpp"

//#include <cmath> // for fabs(x)

/******************************************************************************/

Boundary::Boundary(PLBPD *plbpd): Pointers(plbpd){
  locforce = NULL;
  loctorque = NULL;
  buffer = NULL;
  xperiodic = false;
  yperiodic = false;
  zperiodic = false;
}

/******************************************************************************/

Boundary::~Boundary(){
  if(locforce) delete [] locforce;
  if(loctorque) delete [] loctorque;
  if(buffer) delete [] buffer;
}

/******************************************************************************/

void Boundary::setup(){

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  // called 'after' domain->setup() and lattice->setup()
  for(int i=0; i<6; i++){
    type[i] = input->bdryType[i];
  }
  if(type[0]=="p") xperiodic = true;
  if(type[2]=="p") yperiodic = true;
  if(type[4]=="p") zperiodic = true;

  for(int i=0; i<3; i++){
    uE[i] = input->uE[i];
    uW[i] = input->uW[i];
    uN[i] = input->uN[i];
    uS[i] = input->uS[i];
    uU[i] = input->uU[i];
    uD[i] = input->uD[i];
  }
  
  /*
    set boundary flags
    in our simulations, nodes on the west boundary (not including edges)
    can be 'inflow (fluid)' or obstacle (noslip/freeslip). nodes on
    the east boundary (not including edges) can be 'outflow (fluid)' or
    obstacle (noslip).
    the north/south boundaries (including edges in the x directions,
    but not along the z direction) and up/down boundaries (includes edges
    in the x and y directions) are obstacle nodes (noslip)
  */

  int intfld, inflow, outflow, noslip;
  intfld = inflow = outflow = noslip = 0;

  // for domain nodes only
  for(int zl=2; zl<zdim+2; zl++){
    for(int yl=2; yl<ydim+2; yl++){
      for(int xl=2; xl<xdim+2; xl++){
	// global coords, 1<=xg<=Lx etc.
	int xg = (xl - 1) + domain->lxmin;
	int yg = (yl - 1) + domain->lymin;
	int zg = (zl - 1) + domain->lzmin;

	// east boundary, can be outflow, noslip, periodic
	if( (xg==domain->Lx) and 
	    ((yg>=1) and (yg<=domain->Ly)) and 
	    ((zg>=1) and (zg<=domain->Lz)) and
	    (type[0]!="p") ){
	  for(int i=0; i<3; i++)
	    lattice->u[zl][yl][xl][i] = uE[i];
	  if(type[0]=="o"){
	    lattice->type[zl][yl][xl] = OUTFLOW;
	    outflow++;
	  }
	  else if(type[0]=="n"){
	    lattice->type[zl][yl][xl] = NOSLIP;
	    noslip++;
	  }
	  else{
	    string str("invalid EAST boundary type");
	    error->all(str);
	  }
	}
	  
	// west boundary, can be inflow, noslip, periodic
	if( (xg==1) and 
	    ((yg>=1) and (yg<=domain->Ly)) and 
	    ((zg>=1) and (zg<=domain->Lz)) and
	    (type[1]!="p") ){
	  for(int i=0; i<3; i++)
	    lattice->u[zl][yl][xl][i] = uW[i];
	  if(type[1]=="i"){
	    lattice->type[zl][yl][xl] = INFLOW;
	    inflow++;
	  }
	  else if(type[1]=="n"){
	    lattice->type[zl][yl][xl] = NOSLIP;
	    noslip++;
	  }
	  else{
	    string str("invalid WEST boundary type");
	    error->all(str);
	  }
	}

	// north boundary, can be noslip, periodic
	if( (yg==domain->Ly) and
	    ((xg>=1) and (xg<=domain->Lx)) and 
	    ((zg>=1) and (zg<=domain->Lz)) and
	    (type[2]!="p") ){
	  for(int i=0; i<3; i++)
	    lattice->u[zl][yl][xl][i] = uN[i];
	  if(type[2]=="n"){
	    lattice->type[zl][yl][xl] = NOSLIP;
	    noslip++;
	  }
	  else{
	    string str("invalid NORTH boundary type");
	    error->all(str);
	  }
	}
	
	// south boundary, can be noslip, periodic
	if( (yg==1) and
	    ((xg>=1) and (xg<=domain->Lx)) and 
	    ((zg>=1) and (zg<=domain->Lz)) and
	    (type[3]!="p") ){
	  for(int i=0; i<3; i++)
	    lattice->u[zl][yl][xl][i] = uS[i];
	  if(type[3]=="n"){
	    lattice->type[zl][yl][xl] = NOSLIP;
	    noslip++;
	  }
	  else{
	    string str("invalid SOUTH boundary type");
	    error->all(str);
	  }
	}

	// up boundary, can be noslip, periodic
	if( (zg==domain->Lz) and
	    ((yg>=1) and (yg<=domain->Ly)) and 
	    ((xg>=1) and (xg<=domain->Lx)) and
	    (type[4]!="p") ){
	  for(int i=0; i<3; i++)
	    lattice->u[zl][yl][xl][i] = uU[i];
	  if(type[4]=="n"){
	    lattice->type[zl][yl][xl] = NOSLIP;
	    noslip++;
	  }
	  else{
	    string str("invalid UP boundary type");
	    error->all(str);
	  }
	}
	
	// down boundary, can be freeslip or noslip
	if( (zg==1) and
	    ((yg>=1) and (yg<=domain->Ly)) and 
	    ((xg>=1) and (xg<=domain->Lx)) and
	    (type[5]!="p") ){
	  for(int i=0; i<3; i++)
	    lattice->u[zl][yl][xl][i] = uD[i];
	  if(type[5]=="n"){
	    lattice->type[zl][yl][xl] = NOSLIP;
	    noslip++;
	  }
	  else{
	    string str("invalid DOWN boundary type");
	    error->all(str);
	  }
	}
      }
    }
  }

  // allocate memory for locforce, loctorque and buffer
  // to compute the hydrodynamic force on the agglomerate(s)
  locforce = new double[3*input->nAg];
  loctorque = new double[3*input->nAg];
  buffer = new double[3*input->nAg];
}

/******************************************************************************/

void Boundary::setup_heatxfer(){
  // called 'after' boundary->setup()
  // local copy of (sub) domain dimensions
  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  /*
    in case of thermal bdry nodes, the corners are treated as follows:
    suppose we have a 2d case, E/W are dirichlet bdries and N/S are neumann.
    then, both dirichletW and neumannS would contain the bottom-left node, which
    by design is a neumann node. now since the neumann conditions are applied
    after the dirichlet condns, the bottom-left node would satify the
    neumann condition
  */

  // for domain and communicated nodes 
  for(int zl=1; zl<zdim+3; zl++){
    for(int yl=1; yl<ydim+3; yl++){
      for(int xl=1; xl<xdim+3; xl++){
	// global coords, 1<=xg<=Lx etc.
	int xg = (xl - 1) + domain->lxmin;
	int yg = (yl - 1) + domain->lymin;
	int zg = (zl - 1) + domain->lzmin;

	if( (xg==domain->Lx) and 
	    ((yg>0) and (yg<domain->Ly+1)) and 
	    ((zg>0) and (zg<domain->Lz+1)) and
	    (type[0]=="n") ){
	  if(input->bdryType_th[0]=="d"){ // dirichlet east
	    tuples3<int> node(xl,yl,zl);
	    dirichlet[0].push_back(node);
	  }
	  else if(input->bdryType_th[0]=="n"){
	    tuples3<int> node(xl,yl,zl);
	    neumann[0].push_back(node);
	  }
	  else{
	    string str = "invalid E thermal bdry condn";
	    error->all(str);
	  }
	}
	
	if( (xg==1) and 
	    ((yg>0) and (yg<domain->Ly+1)) and 
	    ((zg>0) and (zg<domain->Lz+1)) and
	    (type[1]=="n") ){
	  if(input->bdryType_th[1]=="d"){ // dirichlet west
	    tuples3<int> node(xl,yl,zl);
	    dirichlet[1].push_back(node);
	  }
	  else if(input->bdryType_th[1]=="n"){
	    tuples3<int> node(xl,yl,zl);
	    neumann[1].push_back(node);
	  }
	  else{
	    string str = "invalid W thermal bdry condn";
	    error->all(str);
	  }
	}

	if( (yg==domain->Ly) and
	    ((xg>0) and (xg<domain->Lx+1)) and 
	    ((zg>0) and (zg<domain->Lz+1)) and
	    (type[2]=="n") ){
	  if(input->bdryType_th[2]=="d"){ // dirichlet north
	    tuples3<int> node(xl,yl,zl);
	    dirichlet[2].push_back(node);
	  }
	  else if(input->bdryType_th[2]=="n"){
	    tuples3<int> node(xl,yl,zl);
	    neumann[2].push_back(node);
	  }
	  else{
	    string str = "invalid N thermal bdry condn";
	    error->all(str);
	  }
	}

	// south boundary: if noslip can be either dirichlet or neumann
	if( (yg==1) and
	    ((xg>0) and (xg<domain->Lx+1)) and 
	    ((zg>0) and (zg<domain->Lz+1)) and
	    (type[3]=="n") ){
	  if(input->bdryType_th[3]=="d"){ // dirichlet south
	    tuples3<int> node(xl,yl,zl);
	    dirichlet[3].push_back(node);
	  }
	  else if(input->bdryType_th[3]=="n"){
	    tuples3<int> node(xl,yl,zl);
	    neumann[3].push_back(node);
	  }
	  else{
	    string str = "invalid S thermal bdry condn";
	    error->all(str);
	  }
	}

	// up boundary: if noslip, can be either dirichlet or neumann
	if( (zg==domain->Lz) and
	    ((yg>0) and (yg<domain->Ly+1)) and 
	    ((xg>0) and (xg<domain->Lx+1)) and
	    (type[4]=="n") ){
	  if(input->bdryType_th[4]=="d"){ // dirichlet up
	    tuples3<int> node(xl,yl,zl);
	    dirichlet[4].push_back(node);
	  }
	  else if(input->bdryType_th[4]=="n"){
	    tuples3<int> node(xl,yl,zl);
	    neumann[4].push_back(node);
	  }
	  else{
	    string str = "invalid U thermal bdry condn";
	    error->all(str);
	  }
	} 

	// down boundary: if noslip, can be either dirichlet or neumann
	if( (zg==1) and
	    ((yg>0) and (yg<domain->Ly+1)) and 
	    ((xg>0) and (xg<domain->Lx+1)) and
	    (type[5]=="n") ){
	  if(input->bdryType_th[5]=="d"){ // dirichlet down
	    tuples3<int> node(xl,yl,zl);
	    dirichlet[5].push_back(node);
	  }
	  else if(input->bdryType_th[5]=="n"){
	    tuples3<int> node(xl,yl,zl);
	    neumann[5].push_back(node);
	  }
	  else{
	    string str = "invalid D thermal bdry condn";
	    error->all(str);
	  }
	}
      }
    }
  }

  /*
  for(int direction=0; direction<6; direction++){
    cout<<"dirn: "<<direction<<endl;
    cout<<"dirichlet\n";
    for(int iNode=0; iNode<dirichlet[direction].size(); iNode++)
      cout<<dirichlet[direction][iNode]<<endl;
    cout<<"neumann\n";
    for(int iNode=0; iNode<neumann[direction].size(); iNode++)
      cout<<neumann[direction][iNode]<<endl;
  }
  */
}

/******************************************************************************/

void Boundary::apply_fluidBC(){
  timer->start(APPLYBC_FL);

  int xdim = domain->lx;
  int ydim = domain->ly;
  int zdim = domain->lz;

  int nx, ny, nz;    // nbring nodes

  int nAg = p_dynamics->nAglmrts;
  Agglomerate *aglmrt = p_dynamics->aglmrt;
  double ****nloc = lattice->n;
  double ****u = lattice->u;
  double ***rho = lattice->rho;
  int ***nodetype = lattice->type;

  double cu, rb[3], ub[3], force[3];

  int nV = model->fluid->nVelocity;
  double cs2 = model->fluid->cs2;
  int (*c)[3] = model->fluid->c;
  double *w = model->fluid->w;
  int *reverse = model->fluid->reverse;

  // first initialize all (local) forces and torques on agglomerate(s)
  for(int iAg=0; iAg<nAg; iAg++){
    locforce[3*iAg+0]  = 0.;
    locforce[3*iAg+1]  = 0.;
    locforce[3*iAg+2]  = 0.;
    loctorque[3*iAg+0] = 0.;
    loctorque[3*iAg+1] = 0.;
    loctorque[3*iAg+2] = 0.;
  }    
  
  // domain and bdry layer nodes only
  for(int zl=1; zl<zdim+3; zl++){
    for(int yl=1; yl<ydim+3; yl++){
      for(int xl=1; xl<xdim+3; xl++){
	// NOSLIP node
	if(nodetype[zl][yl][xl]==NOSLIP){

	  // fluid boundary condition
	  for(int i=0; i<3; i++)
	    ub[i] = u[zl][yl][xl][i];       // bdry velocity

	  for(int kp=1; kp<nV; kp++){       // kp=0 gives this node coords
	    nx = xl + c[kp][0];             // nbring coords
	    ny = yl + c[kp][1];
	    nz = zl + c[kp][2];
	    
	    if((nodetype[nz][ny][nx]==INTFLD) or
	       (nodetype[nz][ny][nx]==INFLOW) or
	       (nodetype[nz][ny][nx]==OUTFLOW)){
	      // nbr is fluid
	      int k = reverse[kp];
	      cu = c[k][0]*ub[0] + c[k][1]*ub[1] + c[k][2]*ub[2];
	      nloc[nz][ny][nx][kp] = nloc[zl][yl][xl][k] - 
		2.*w[k]*rho[nz][ny][nx]*cu/cs2;
	    }
	    else if(nodetype[nz][ny][nx]>=0){
	      // nbr is aglmrt solid
	      string str("aglmrt too close to NOSLIP node");
	      error->all(str);
	    }
	  }
	}

	// AGLMRT node
	else if(nodetype[zl][yl][xl]>=0){   // aglmrt types are >=0

	  int iAg = nodetype[zl][yl][xl];   // aglmrt ID
	  for(int kp=1; kp<nV; kp++){
	    nx = xl + c[kp][0];
	    ny = yl + c[kp][1];
	    nz = zl + c[kp][2];
	    int jAg = nodetype[nz][ny][nx];
	    
	    if(jAg==INTFLD){                
	      // nbring node in interior fluid
	      int k = reverse[kp];
	      calc_ub(xl,yl,zl,kp,rb,ub,aglmrt,iAg);
	      cu = c[k][0]*ub[0] + c[k][1]*ub[1] + c[k][2]*ub[2];
	      // update dist at the fluid node
	      double tmp = w[k]*rho[nz][ny][nx]*cu/cs2;
	      nloc[nz][ny][nx][kp] = nloc[zl][yl][xl][k] - 2.*tmp;

	      // compute force and torque
	      if( (xl>1) and (xl<xdim+2) and
	      	  (yl>1) and (yl<ydim+2) and
	      	  (zl>1) and (zl<zdim+2) ){
		// do NOT to count the force on the boundary nodes twice).
		// force
		force[0] = 2.*(nloc[nz][ny][nx][kp] + tmp)*c[k][0];
		force[1] = 2.*(nloc[nz][ny][nx][kp] + tmp)*c[k][1];
		force[2] = 2.*(nloc[nz][ny][nx][kp] + tmp)*c[k][2];

		// NOTE: the moment particles are introduced, small
		// errors creep up between serial and parallel versions
		// because of floating point rounding errors
		// e.g -0.1111111111111111327 + -0.2222222222222222654
		// = -0.3333333333333333703 instead of -0.3333333333333333981.
		// in case of parallel processing if this sum is avoided
		// in one of the procs, the error does not rise.
		locforce[3*iAg+0] += force[0];
		locforce[3*iAg+1] += force[1];
		locforce[3*iAg+2] += force[2];

		/*
		// debug
		int xg = (xl - 2) + domain->lxmin;
		int yg = (yl - 2) + domain->lymin;
		int zg = (zl - 2) + domain->lzmin;
		int nxg = (nx - 2) + domain->lxmin;
		int nyg = (ny - 2) + domain->lymin;
		int nzg = (nz - 2) + domain->lzmin;
		printf("(%2d,%2d,%2d):(%2d,%2d,%2d), %22.19f %22.19f %22.19f\n",
		       xg, yg, zg, nxg, nyg, nzg,
		       force[3*iAg+0], force[3*iAg+1], force[3*iAg+2]);
		printf("%22.19f %22.19f %22.19f\n", locforce[3*iAg+0],
		       locforce[3*iAg+1], locforce[3*iAg+2]);
		*/

		// torque = rb x force
		// loctorque[3*iAg+0] += rb[1]*force[2] - rb[2]*force[1];
		// loctorque[3*iAg+1] += rb[2]*force[0] - rb[0]*force[2];
		// loctorque[3*iAg+2] += rb[0]*force[1] - rb[1]*force[0];

		double torque[] = {0.,0.,0.};
		torque[0] += rb[1]*force[2] - rb[2]*force[1];
		torque[1] += rb[2]*force[0] - rb[0]*force[2];
		torque[2] += rb[0]*force[1] - rb[1]*force[0];

		loctorque[3*iAg+0] += torque[0];
		loctorque[3*iAg+1] += torque[1];
		loctorque[3*iAg+2] += torque[2];

		// printf("(%d,%d,%d): %d: (%e,%e,%e), (%e,%e,%e), (%e,%e,%e)\n", 
		//        xl,yl,zl,kp,rb[0],rb[1],rb[2],force[0],force[1],force[2],
		//        torque[0],torque[1],torque[2]);
	      }
	    }
	    
	    else if((jAg>=0) and (jAg!=iAg)){
	      string str = "aglmrts " + utils::T2str<int>(iAg) + " and " +
		utils::T2str<int>(jAg) + " are too close";
	      str += " (" + utils::T2str<int>(xl) + "," 
		+ utils::T2str<int>(yl) + "," 
		+ utils::T2str<int>(zl) + "), (" 
		+ utils::T2str<int>(nx) + "," 
		+ utils::T2str<int>(ny) + "," 
		+ utils::T2str<int>(nz) + ")";
	      error->all(str);
	    }

	    else if((jAg==INFLOW) or (jAg==OUTFLOW)){
	      string str = "aglmrt " + utils::T2str<int>(iAg) + " too close " +
		"to INFLOW/OUTFLOW bdry";
	      error->all(str);
	    }
	  }
	}
	
	// inflow (some calculations are repeated for each inflow node)
	else if(nodetype[zl][yl][xl]==INFLOW){  // ONLY at the WEST boundary
	  lb_dynamics->calc_eqdist(rho[zl][yl][xl], uW, nloc[zl][yl][xl]);
	}

	// outflow
	else if(nodetype[zl][yl][xl]==OUTFLOW){ // ONLY at the EAST boundary
	  lb_dynamics->calc_eqdist(rho[zl][yl][xl], uE, nloc[zl][yl][xl]);
	}
      } // end for
    } // end for
  } // end for

  if(nAg>0){
    // add force across procs
    // i need force on all procs (hence Allreduce)
    int count = 3*nAg;
    MPI_Allreduce(locforce,buffer,count,MPI_DOUBLE,MPI_SUM,world);
    // copy force (average of the current and the previous values)
    // to p_dynamics->aglmrt[iAg].F_hyd
    for(int iAg=0; iAg<nAg; iAg++){
      aglmrt[iAg].F_hyd[0] = 0.5*(buffer[3*iAg+0] + aglmrt[iAg].F_hyd[0]);
      aglmrt[iAg].F_hyd[1] = 0.5*(buffer[3*iAg+1] + aglmrt[iAg].F_hyd[1]);
      aglmrt[iAg].F_hyd[2] = 0.5*(buffer[3*iAg+2] + aglmrt[iAg].F_hyd[2]);
      // aglmrt[iAg].F_hyd[0] = buffer[3*iAg+0];
      // aglmrt[iAg].F_hyd[1] = buffer[3*iAg+1];
      // aglmrt[iAg].F_hyd[2] = buffer[3*iAg+2];
    }

    MPI_Barrier(world);

    // add torque across procs
    MPI_Allreduce(loctorque,buffer,count,MPI_DOUBLE,MPI_SUM,world);
    // copy torque (average of the current and the previous values)
    // to p_dynamics->aglmrt[iAg].T_hyd
    for(int iAg=0; iAg<nAg; iAg++){
      aglmrt[iAg].T_hyd[0] = 0.5*(buffer[3*iAg+0] + aglmrt[iAg].T_hyd[0]);
      aglmrt[iAg].T_hyd[1] = 0.5*(buffer[3*iAg+1] + aglmrt[iAg].T_hyd[1]);
      aglmrt[iAg].T_hyd[2] = 0.5*(buffer[3*iAg+2] + aglmrt[iAg].T_hyd[2]);
      // aglmrt[iAg].T_hyd[0] = buffer[3*iAg+0];
      // aglmrt[iAg].T_hyd[1] = buffer[3*iAg+1];
      // aglmrt[iAg].T_hyd[2] = buffer[3*iAg+2];
      // if(comm->me==0){
      // 	cout<<iAg<<": ("<<aglmrt[iAg].F_hyd[0]<<","<<aglmrt[iAg].F_hyd[1]
      // 	    <<","<<aglmrt[iAg].F_hyd[2]<<"), (";
      // 	cout<<aglmrt[iAg].T_hyd[0]<<","<<aglmrt[iAg].T_hyd[1]
      // 	    <<","<<aglmrt[iAg].T_hyd[2]<<")\n";      
      // } 
    }
  }
  
  timer->end(APPLYBC_FL);
}

/******************************************************************************/

void Boundary::calc_ub(int xs, int ys, int zs, int k, double *rb,
		       double *ub, Agglomerate *aglmrt, int iAg){

  /* calculate \vec{u_b}, the velocity of the solid (aglmrt) iAg at
     the boundary (xs,ys,zs)+0.5*\vec{c_K}. \vec{c_k} is the direction
     towards the boundary from the solid node (xs,ys,zs) so that the 
     boundary node is halfway between (xs,ys,zs) and the node in the
     direction \vec{c_k} 
     \vec{r_b} is the vector pointing from the com of the aglmrt iAg
     to the boudary node
  */
  
  double Lx = static_cast<double>(domain->Lx);
  double Ly = static_cast<double>(domain->Ly);
  double Lz = static_cast<double>(domain->Lz);  

  // Reacll: x=2, y=2, z=2 in proc 0 has coords (0,0,0)
  double xsg, ysg, zsg;  // box coordinates
  xsg = static_cast<double>((xs - 2) + domain->lxmin);
  ysg = static_cast<double>((ys - 2) + domain->lymin);
  zsg = static_cast<double>((zs - 2) + domain->lzmin);

  // boundary node coordinates (bdry is in the direction k)
  int (*c)[3] = model->fluid->c;
  double tmp[] = {0., 0., 0.};
  tmp[0] = xsg + 0.5*static_cast<double>(c[k][0]);
  tmp[1] = ysg + 0.5*static_cast<double>(c[k][1]);
  tmp[2] = zsg + 0.5*static_cast<double>(c[k][2]);
  // this is needed to put coords like tmp[0] = -0.5 back in the box, 39.5
  if(xperiodic) tmp[0] = utils::mod(tmp[0],Lx);
  if(yperiodic) tmp[1] = utils::mod(tmp[1],Ly);
  if(zperiodic) tmp[2] = utils::mod(tmp[2],Lz);
  
  /* boundary node velocity
     ub = U + Omega x (tmp-X) .. Eqn 2.8, Ladd 1994b
     where:
       U:     lin vel of the com of the aglmrt iAg
       Omega: ang vel about the com of iAg
       X:     com of iAg */
  double X[3], *U=NULL, *Omega=NULL;
  X[0] = aglmrt[iAg].com[0];
  X[1] = aglmrt[iAg].com[1];
  X[2] = aglmrt[iAg].com[2];
  // correction for periodic bdry
  if(xperiodic) X[0] = utils::mod(X[0],Lx);  
  if(yperiodic) X[1] = utils::mod(X[1],Ly);
  if(zperiodic) X[2] = utils::mod(X[2],Lz);

  U = aglmrt[iAg].u;
  Omega = aglmrt[iAg].omega;

  // rb = tmp - X
  for(int i=0; i<3; i++)
    rb[i] = tmp[i] - X[i];

  // correction for periodic bdry
  if(xperiodic) rb[0] -= utils::nint(rb[0]/Lx)*Lx;
  if(yperiodic) rb[1] -= utils::nint(rb[1]/Ly)*Ly;
  if(zperiodic) rb[2] -= utils::nint(rb[2]/Lz)*Lz;

  // compute ub
  ub[0] = U[0] + Omega[1]*rb[2] - Omega[2]*rb[1];
  ub[1] = U[1] - Omega[0]*rb[2] + Omega[2]*rb[0];
  ub[2] = U[2] + Omega[0]*rb[1] - Omega[1]*rb[0];

  // printf("(%19.16e,%19.16e,%19.16e):(%19.16e,%19.16e,%19.16e): %d  (%e,%e,%e)\n", 
  // 	 tmp[0],tmp[1],tmp[2],X[0],X[1],X[2],k,rb[0],rb[1],rb[2]);
  
}

/******************************************************************************/

void Boundary::apply_thermalBC(){
  timer->start(APPLYBC_TH);

  if(!input->heat_xfer){
    string str = "Boundary::apply_thermalBC() called even though";
    str += " heat transfer is not part of the simulation";
    error->all(str);
  }

  double *** const theta = lattice->theta;
  double **** const nth  = lattice->nth;
  int nV_th = model->thermal->nVelocity;

  // first all dirichlet nodes (in order E, W, N, S, U, D)
  double btheta[6], *uBdry = NULL;
  for(int direction=0; direction<6; direction++){
    btheta[direction] = input->bdryTemp[direction];
    if     (direction==0) uBdry = uE;
    else if(direction==1) uBdry = uW;
    else if(direction==2) uBdry = uN;
    else if(direction==3) uBdry = uS;
    else if(direction==4) uBdry = uU;
    else if(direction==5) uBdry = uD;
    else{
      string str = "invalid value of 'direction'";
      error->all(str);
    }
    for(int iNode=0; iNode<dirichlet[direction].size(); iNode++){
      tuples3<int> node = dirichlet[direction][iNode];
      int xl = node.a; int yl = node.b; int zl = node.c;
      th_dynamics->calc_eqdist(btheta[direction],uBdry,nth[zl+1][yl+1][xl+1]);
    }
  }
    
  // next all neumann nodes (in order as above)
  double bdtheta[6];
  for(int direction=0; direction<6; direction++){
    bdtheta[direction] = input->bdryTemp[direction];
    for(int iNode=0; iNode<neumann[direction].size(); iNode++){
      tuples3<int> node = neumann[direction][iNode];
      int xl = node.a; int yl = node.b; int zl = node.c;
      if(direction==0){
	for(int k=0; k<nV_th; k++)
	  nth[zl+1][yl+1][xl+1][k] = nth[zl+1][yl+1][xl+1-1][k];
      }
      if(direction==1){
	for(int k=0; k<nV_th; k++)
	  nth[zl+1][yl+1][xl+1][k] = nth[zl+1][yl+1][xl+1+1][k];
      }
      if(direction==2){
	for(int k=0; k<nV_th; k++)
	  nth[zl+1][yl+1][xl+1][k] = nth[zl+1][yl+1-1][xl+1][k];
      }
      if(direction==3){
	for(int k=0; k<nV_th; k++)
	  nth[zl+1][yl+1][xl+1][k] = nth[zl+1][yl+1+1][xl+1][k];
      }
      if(direction==4){
	for(int k=0; k<nV_th; k++)
	  nth[zl+1][yl+1][xl+1][k] = nth[zl+1-1][yl+1][xl+1][k];
      }
      if(direction==5){
	for(int k=0; k<nV_th; k++)
	  nth[zl+1][yl+1][xl+1][k] = nth[zl+1+1][yl+1][xl+1][k];
      }
    }
  }
  timer->end(APPLYBC_TH);
}

/******************************************************************************/
