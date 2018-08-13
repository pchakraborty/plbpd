#include "SimData.hpp"

#include <mm_malloc.h>
#include <string>
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

SimData::SimData(const std::tuple<size_t, size_t, size_t> domain_dimensions,
                 const size_t num_directions)
    : _domain_dimensions(domain_dimensions), _kdim(num_directions) {
    rho = nullptr;
    u = nullptr;
    n = nullptr;
    ntmp = nullptr;
    _bootstrap();
}

SimData::~SimData() {}

void SimData::_bootstrap() {
    /*
      The variables n and ntmp require a buffer layer on each side
      (for streaming), so size(n) = size(ntmp) = (_zdim+2, _ydim+2, _xdim+2)
      Variables, nodetype, rho, u can very well have dimensions
      (_zdim, _ydim, _xdim). However, for uniformity of indexing, all the
      arrays have the same dimensions.

      0:_xdim+1 include the buffer layer, where
      1:_xdim   does not include the buffer layer
    */
    size_t _xdim, _ydim, _zdim;
    std::tie(_xdim, _ydim, _zdim) = _domain_dimensions;

    // rho = new array3f(_zdim+2, _ydim+2, _xdim+2, 0.0f);
    // u = new array4f(_zdim+2, _ydim+2, _xdim+2, 3, 0.0f);
    // n = new array4f(_zdim+2, _ydim+2, _xdim+2, _kdim, 0.0f);
    // ntmp = new array4f(_zdim+2, _ydim+2, _xdim+2, _kdim, 0.0f);

    rho = new Field::ScalarField<float, 1>(_zdim, _ydim, _xdim, 0.0f);
    u = new Field::VectorField<float, 1>(_zdim, _ydim, _xdim, 3, 0.0f);
    n = new Field::VectorField<float, 1>(_zdim, _ydim, _xdim, _kdim, 0.0f);
    ntmp = new Field::VectorField<float, 1>(_zdim, _ydim, _xdim, _kdim, 0.0f);

    // // allocate aligned mem via _mm_malloc
    // const auto alignment = 64;

    // auto rhosize = (_zdim+2)*(_ydim+2)*(_xdim+2)*sizeof(float);
    // assert(rhosize%alignment==0);
    // rho = static_cast<float*>(_mm_malloc(rhosize, alignment));
    // for (auto i=0; i<(_zdim+2)*(_ydim+2)*(_xdim+2); ++i)
    //     rho[i] = 1.0f;

    // auto usize = (_zdim+2)*(_ydim+2)*(_xdim+2)*3*sizeof(float);
    // assert(usize%alignment==0);
    // u = static_cast<float*>(_mm_malloc(usize, alignment));
    // for (auto i=0; i<(_zdim+2)*(_ydim+2)*(_xdim+2)*3; ++i)
    //     u[i] = 0.0f;

    // auto nsize = (_zdim+2)*(_ydim+2)*(_xdim+2)*_kdim*sizeof(float);
    // assert(nsize%alignment==0);
    // n = static_cast<float*>(_mm_malloc(nsize, alignment));
    // ntmp = static_cast<float*>(_mm_malloc(nsize, alignment));
    // for (auto i=0; i<(_zdim+2)*(_ydim+2)*(_xdim+2)*_kdim; ++i){
    //     n[i] = 0.0f;
    //     ntmp[i] = 0.0f;
    // }
}

void SimData::write_state(std::string dump_file) {
    auto file = H5Fcreate(dump_file.c_str(),
                          H5F_ACC_TRUNC,
                          H5P_DEFAULT,
                          H5P_DEFAULT);

    // Domain dimensions
    size_t _xdim, _ydim, _zdim;
    std::tie(_xdim, _ydim, _zdim) = _domain_dimensions;

    // Data spaces
    std::array<hsize_t, 3> dims_rho = {_zdim+2, _ydim+2, _xdim+2};
    std::array<hsize_t, 4> dims_u = {_zdim+2, _ydim+2, _xdim+2, 3};
    // std::array<hsize_t, 4> dims_n = {_zdim+2, _ydim+2, _xdim+2, _kdim};

    auto dataspace_rho = H5Screate_simple(3, dims_rho.data(), NULL);
    auto dataspace_u = H5Screate_simple(4, dims_u.data(), NULL);
    // auto dataspace_n = H5Screate_simple(4, dims_n.data(), NULL);

    // Data set
    auto dataset_rho = H5Dcreate2(
        file,
        "/Density",
        H5T_NATIVE_FLOAT,
        dataspace_rho,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT);
    auto dataset_u = H5Dcreate2(
        file,
        "/FlowVelocity",
        H5T_NATIVE_FLOAT,
        dataspace_u,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT);
    // auto dataset_n = H5Dcreate2(
    //     file,
    //     "/ParticleDistribution",
    //     H5T_NATIVE_FLOAT,
    //     dataspace_n,
    //     H5P_DEFAULT,
    //     H5P_DEFAULT,
    //     H5P_DEFAULT);

    herr_t status;

    // Write data set
    status = H5Dwrite(
        dataset_rho,
        H5T_NATIVE_FLOAT,
        H5S_ALL,
        H5S_ALL,
        H5P_DEFAULT,
        rho->get(0, 0, 0));
    status = H5Dwrite(
        dataset_u,
        H5T_NATIVE_FLOAT,
        H5S_ALL, H5S_ALL,
        H5P_DEFAULT,
        u->get(0, 0, 0, 0));
    // status = H5Dwrite(
    //     dataset_n,
    //     H5T_NATIVE_FLOAT,
    //     H5S_ALL,
    //     H5S_ALL,
    //     H5P_DEFAULT,
    //     n->get(0, 0, 0, 0));

    // clean up
    status = H5Dclose(dataset_rho);
    status = H5Dclose(dataset_u);
    // status = H5Dclose(dataset_n);
    status = H5Sclose(dataspace_rho);
    status = H5Sclose(dataspace_u);
    // status = H5Sclose(dataspace_n);
    status = H5Fclose(file);
  }

void SimData::_read_state() {
    throw std::logic_error("_read_state() has not yet been implemented");
}
