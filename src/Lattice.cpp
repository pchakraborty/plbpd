#include "Lattice.hpp"
#include <mm_malloc.h>
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
    const auto kdim = _lbmodel->getNumberOfDirections();
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim)  = _domain->getDimensions();
    
    rho = new array3f(zdim+2, ydim+2, xdim+2);
    u = new array4f(zdim+2, ydim+2, xdim+2, 3);
    n = new array4f(zdim+2, ydim+2, xdim+2, kdim);
    ntmp = new array4f(zdim+2, ydim+2, xdim+2, kdim);
    
    // // allocate mem via new
    // rho = new float[(zdim+2)*(ydim+2)*(xdim+2)](); // zero initialized via ()
    // u = new float[(zdim+2)*(ydim+2)*(xdim+2)*3]();
    // n = new float[(zdim+2)*(ydim+2)*(xdim+2)*kdim]();
    // ntmp = new float[(zdim+2)*(ydim+2)*(xdim+2)*kdim]();

    // // allocate aligned mem via _mm_malloc
    // const auto alignment = 64;

    // auto rhosize = (zdim+2)*(ydim+2)*(xdim+2)*sizeof(float);
    // assert(rhosize%alignment==0);
    // rho = static_cast<float*>(_mm_malloc(rhosize, alignment));
    // for (auto i=0; i<(zdim+2)*(ydim+2)*(xdim+2); ++i)
    //     rho[i] = 1.0f;
    
    // auto usize = (zdim+2)*(ydim+2)*(xdim+2)*3*sizeof(float);
    // assert(usize%alignment==0);
    // u = static_cast<float*>(_mm_malloc(usize, alignment));
    // for (auto i=0; i<(zdim+2)*(ydim+2)*(xdim+2)*3; ++i)
    //     u[i] = 0.0f;
    
    // auto nsize = (zdim+2)*(ydim+2)*(xdim+2)*kdim*sizeof(float);
    // assert(nsize%alignment==0);
    // n = static_cast<float*>(_mm_malloc(nsize, alignment));
    // ntmp = static_cast<float*>(_mm_malloc(nsize, alignment));
    // for (auto i=0; i<(zdim+2)*(ydim+2)*(xdim+2)*kdim; ++i){
    //     n[i] = 0.0f;
    //     ntmp[i] = 0.0f;
    // }
}

void Lattice::writeState(std::string dumpFile){
    auto file = H5Fcreate(dumpFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Data spaces
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim)  = _domain->getDimensions();
    auto kdim = _lbmodel->getNumberOfDirections();

    std::array<hsize_t, 3> dims_rho = {zdim+2, ydim+2, xdim+2};
    std::array<hsize_t, 4> dims_u = {zdim+2, ydim+2, xdim+2, 3};
    // std::array<hsize_t, 4> dims_n = {zdim+2, ydim+2, xdim+2, kdim};

    auto dataspace_rho = H5Screate_simple(3, dims_rho.data(), NULL);
    auto dataspace_u = H5Screate_simple(4, dims_u.data(), NULL);
    // auto dataspace_n = H5Screate_simple(4, dims_n.data(), NULL);
    
    // Data set
    auto dataset_rho = H5Dcreate2(file, "/Density", H5T_NATIVE_FLOAT, dataspace_rho, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto dataset_u = H5Dcreate2(file, "/FlowVelocity", H5T_NATIVE_FLOAT, dataspace_u, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // auto dataset_n = H5Dcreate2(file, "/ParticleDistribution", H5T_NATIVE_FLOAT, dataspace_n, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    herr_t status;

    // Write data set
    status = H5Dwrite(dataset_rho, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho->get());
    status = H5Dwrite(dataset_u, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, u->get());
    // status = H5Dwrite(dataset_n, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n->get());
    
    // clean up
    status = H5Dclose(dataset_rho);
    status = H5Dclose(dataset_u);
    // status = H5Dclose(dataset_n);
    status = H5Sclose(dataspace_rho);
    status = H5Sclose(dataspace_u);
    // status = H5Sclose(dataspace_n);
    status = H5Fclose(file);
  }

void Lattice::_readState(){
    throw std::logic_error("readState() has not yet been implemented");  
}
