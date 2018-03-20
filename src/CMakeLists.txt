cmake_minimum_required(VERSION 3.10)

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

add_executable(${PROJECT} ${srcFiles})
target_link_libraries(
  ${PROJECT}
  # internal deps (none at this time)
  # external deps
  ${LIBTBB}
  )
