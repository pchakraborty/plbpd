#include "aLattice.hpp"
#include <mm_malloc.h>
#include <iostream>
#include <cassert>
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
    rho = nullptr;
    u = nullptr;
    n = nullptr;
    ntmp = nullptr;
    _bootstrap();
}

Lattice::~Lattice(){
    // // If mem was allocated via new
    // if (rho) delete [] rho;
    // if (rho) delete [] u;
    // if (rho) delete [] n;
    // if (rho) delete [] ntmp;

    // If aligned mem was allocated via _mm_alloc
    if (rho) _mm_free(rho);
    if (u) _mm_free(u);
    if (n) _mm_free(n);
    if (ntmp) _mm_free(ntmp);
}

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
    const auto numVelocityVectors = _lbmodel->getNumVelocityVectors();
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim)  = _domain->getDimensions();

    // // allocate mem via new
    // rho = new float[(zdim+2)*(ydim+2)*(xdim+2)](); // zero initialized via ()
    // u = new float[(zdim+2)*(ydim+2)*(xdim+2)*3]();
    // n = new float[(zdim+2)*(ydim+2)*(xdim+2)*numVelocityVectors]();
    // ntmp = new float[(zdim+2)*(ydim+2)*(xdim+2)*numVelocityVectors]();
    
    // allocate aligned mem via _mm_malloc
    const auto alignment = 64;

    auto rhosize = (zdim+2)*(ydim+2)*(xdim+2)*sizeof(float);
    assert(rhosize%alignment==0);
    rho = static_cast<float*>(_mm_malloc(rhosize, alignment));
    
    auto usize = (zdim+2)*(ydim+2)*(xdim+2)*3*sizeof(float);
    assert(usize%alignment==0);
    u = static_cast<float*>(_mm_malloc(usize, alignment));
    
    auto nsize = (zdim+2)*(ydim+2)*(xdim+2)*numVelocityVectors*sizeof(float);
    assert(nsize%alignment==0);
    n = static_cast<float*>(_mm_malloc(nsize, alignment));
    ntmp = static_cast<float*>(_mm_malloc(nsize, alignment));
}

void Lattice::writeState(){
    auto dumpFile = "state.h5";
    auto file = H5Fcreate(dumpFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Data spaces
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim)  = _domain->getDimensions();
    auto kdim = _lbmodel->getNumVelocityVectors();
    std::array<hsize_t, 4> dims_u = {zdim+2, ydim+2, xdim+2, 3};
    auto dataspace_u = H5Screate_simple(4, dims_u.data(), NULL);
    // std::array<hsize_t, 4> dims_n = {zdim+2, ydim+2, xdim+2, kdim};
    // auto dataspace_n = H5Screate_simple(4, dims_n.data(), NULL);
    
    // Data set
    // auto dataset_n = H5Dcreate2(file, "/ParticleDistribution", H5T_NATIVE_FLOAT,
    //                             dataspace_n, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto dataset_u = H5Dcreate2(file, "/FlowVelocity", H5T_NATIVE_FLOAT,
                                dataspace_u, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    herr_t status;

    // Write data set
    // status = H5Dwrite(dataset_n, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n[0]);
    status = H5Dwrite(dataset_u, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &u[0]);
    
    // clean up
    // status = H5Dclose(dataset_n);
    status = H5Dclose(dataset_u);
    // status = H5Sclose(dataspace_n);
    status = H5Sclose(dataspace_u);
    status = H5Fclose(file);
  }

void Lattice::_readState(){
    throw std::logic_error("readState() has not yet been implemented");  
}
