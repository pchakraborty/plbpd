cmake_minimum_required(VERSION 3.3)

set(exec plbpd)
project(${exec} CXX)
set (CMAKE_CXX_STANDARD 14)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/build-types.cmake)

set(
  srcFiles
  SimData.cpp
  LBModel.cpp
  Domain.cpp
  Boundary.cpp
  Flow.cpp
  Collision.cpp
  CollisionSRT.cpp
  Streaming.cpp
  LBDynamics.cpp
  CalcMoments.cpp
  PLBPD.cpp
  main.cpp
  )

# External dependencies
## TBB
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/tbb.cmake)
## HDF5
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/hdf5.cmake)
include_directories(${HDF5_INCLUDE_DIR})

add_executable(${exec} ${srcFiles})
target_link_libraries(
  ${exec}
  # Dependencies
  ${LIBTBB}
  ${HDF5_C_STATIC_LIBRARY}
  )

# AVX2
option(AVX2 "Build plbpd with AVX2 intrinsics" ON)
if(AVX2)
  message(STATUS "Adding -DAVX2...")
  add_definitions(-DAVX2)
  target_compile_options(${exec} PRIVATE -mfma)
endif()

enable_testing()
add_subdirectory(utest)
