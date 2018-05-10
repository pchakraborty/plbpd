cmake_minimum_required(VERSION 3.2)

set(PROJECT plbpd)
project(${PROJECT} CXX)

set(
  srcFiles
  Lattice.cpp
  LBModel.cpp
  Domain.cpp
  Boundary.cpp
  Flow.cpp
  LBDynamics.cpp
  BGK.cpp
  main.cpp
  )

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
target_compile_options(
  ${PROJECT}
  PRIVATE
  $<$<CONFIG:Release>:-mavx2 -march=core-avx2>
  $<$<CONFIG:Debug>:-mavx2>
  )

enable_testing()
add_subdirectory(utest)
