#include "Lattice.hpp"
#include <iostream>
#include "hdf5.h"

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

Lattice::Lattice(const LBModel *lbmodel, const Domain *domain):
    _lbmodel(lbmodel), _domain(domain){
    _bootstrap();
}

Lattice::~Lattice(){}

void Lattice::_bootstrap(){
    /*
      The variables n and ntmp require a buffer layer on each side
      (for streaming), so size(n) = size(ntmp) = (_zdim+2, _ydim+2, _xdim+2)
      Variables, nodetype, rho, u can very well have dimensions
      (_zdim, _ydim, _xdim). However, for uniformity of indexing, all the arrays
      have the same dimensions.
      
      0:_xdim+1 include the buffer layer, where
      1:_xdim   does not include the buffer layer
    */
    auto numVelocityVectors = _lbmodel->getNumVelocityVectors();
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim)  = _domain->getDimensions();

    nodetype.resize((zdim+2)*(ydim+2)*(xdim+2), 0);
    rho.resize((zdim+2)*(ydim+2)*(xdim+2), 1.0);
    u.resize((zdim+2)*(ydim+2)*(xdim+2)*3, 0.0);
    n.resize((zdim+2)*(ydim+2)*(xdim+2)*numVelocityVectors, 0.0);
    ntmp.resize((zdim+2)*(ydim+2)*(xdim+2)*numVelocityVectors, 0.0);
}

void Lattice::writeState(){
    auto dumpFile = "state.h5";
    auto file = H5Fcreate(dumpFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Data spaces
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim)  = _domain->getDimensions();
    auto kdim = _lbmodel->getNumVelocityVectors();
    std::array<hsize_t, 4> dims_n = {xdim+2, ydim+2, zdim+2, kdim};
    std::array<hsize_t, 4> dims_u = {xdim+2, ydim+2, zdim+2, 3};
    auto dataspace_n = H5Screate_simple(4, dims_n.data(), NULL);
    auto dataspace_u = H5Screate_simple(4, dims_u.data(), NULL);
    
    // Data set
    auto dataset_n = H5Dcreate2(file, "/ParticleDistribution", H5T_NATIVE_FLOAT,
                                dataspace_n, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto dataset_u = H5Dcreate2(file, "/FlowVelocity", H5T_NATIVE_FLOAT,
                                dataspace_u, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    herr_t status;

    // Write data set
    status = H5Dwrite(dataset_n, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n[0]);
    status = H5Dwrite(dataset_u, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &u[0]);
    
    // clean up
    status = H5Dclose(dataset_n);
    status = H5Dclose(dataset_u);
    status = H5Sclose(dataspace_n);
    status = H5Sclose(dataspace_u);
    status = H5Fclose(file);
  }

void Lattice::_readState(){
    throw std::logic_error("readState() has not yet been implemented");  
}
