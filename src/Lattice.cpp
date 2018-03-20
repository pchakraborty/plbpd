#include "Lattice.hpp"
#include <iostream>

/*
  TODO:
  Each MPI task builds it own local domain. In addition to the nodes
  belonging to a given task, the task also contains the boundary nodes
  of the adjacent tasks. Plus, 'n' and 'ntmp' have a buffer layer on
  each side. that makes their dimensions (_zdim+4,_ydim+4,_xdim+4). 
  Variables nodeType, rho, u can very well have dimensions of
  (_zdim+2,_ydim+2,_xdim+2), but that makes the coding difficult. For
  uniformity of indexing, all the arrays have the same dimensions. So,
  0:_xdim+3 includes both the bdry and buffer layers
  1:_xdim+2 includes the bdry layer but excludes the buffer layer
  2:_xdim+1 excludes both the bdry and buffer layers
*/

Lattice::Lattice(const LBModel &lbmodel, const Domain &domain, bool bootstrap=false){
    _numVelocityVectors = lbmodel.getNumVelocityVectors();
    std::tie(_xdim, _ydim, _zdim)  = domain.getDomainDimensions();
    if (bootstrap){
        _bootstrapRestarts();
    } else{
        // TODO: get the name of the restart file
        std::string restartFile = "myRestart.h5";
        _readRestarts(restartFile);
    }
}

Lattice::~Lattice(){}

void Lattice::_bootstrapRestarts(){
    /*
      The variables n and ntmp require a buffer layer on each side
      (for streaming), so size(n) = size(ntmp) = (_zdim+2, _ydim+2, _xdim+2)
      Variables, nodetype, rho, u can very well have dimensions
      For uniformity of indexing, all the arrays have the same dimensions.
      
      0:_xdim+1 include the buffer layer, where
      1:_xdim   does not include the buffer layer
    */
    // allocate and create instances
    nodetype.resize((_zdim+2)*(_ydim+2)*(_xdim+2), 0);
    rho.resize((_zdim+2)*(_ydim+2)*(_xdim+2), 1.0);
    u.resize((_zdim+2)*(_ydim+2)*(_xdim+2)*3, 0.0);
    n.resize((_zdim+2)*(_ydim+2)*(_xdim+2)*_numVelocityVectors, 0.0);
    ntmp.resize((_zdim+2)*(_ydim+2)*(_xdim+2)*_numVelocityVectors, 0.0);
}

void Lattice::_readRestarts(std::string restartFile){
    throw std::logic_error("_readRestarts() has not been implemented");  
}
