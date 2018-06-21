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
  Collision.cpp
  CollisionSRT.cpp
  Streaming.cpp
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
option(AVX2 "Build plbpd with AVX2 intrinsics" ON)
if(AVX2)
  message(STATUS "Adding -DAVX2...")
  add_definitions(-DAVX2)
  target_compile_options(
    ${PROJECT}
    PRIVATE
    $<$<CONFIG:Release>:-mavx2 -march=core-avx2>
    $<$<CONFIG:Debug>:-mavx2>
    )
endif()

enable_testing()
add_subdirectory(utest)
