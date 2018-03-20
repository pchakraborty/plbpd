cmake_minimum_required(VERSION 3.2)

set(PROJECT plbpd)
project(${PROJECT} CXX)

set(srcFiles
  LBModel.cpp
  Domain.cpp
  Lattice.cpp
  LBDynamics.cpp
  BGK.cpp
  main.cpp)

# External dependencies
## TBB
include(./cmake/tbb.cmake)
## HDF5
include(./cmake/hdf5.cmake)
include_directories(${HDF5_INCLUDE_DIR})

add_executable(${PROJECT} ${srcFiles})
target_link_libraries(
  ${PROJECT}
  # Dependencies
  ${LIBTBB}
  ${HDF5_C_STATIC_LIBRARY}
  )
