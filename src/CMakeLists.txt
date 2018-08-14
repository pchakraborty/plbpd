cmake_minimum_required(VERSION 3.3)

set(PROJECT plbpd)
project(${PROJECT} CXX)
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

add_executable(${PROJECT} ${srcFiles})
target_link_libraries(
  ${PROJECT}
  # Dependencies
  ${LIBTBB}
  ${HDF5_C_STATIC_LIBRARY}
  )

# AVX2
option(AVX2 "Build plbpd with AVX2 intrinsics" ON)
if(AVX2)
  message(STATUS "Adding -DAVX2...")
  add_definitions(-DAVX2)
  target_compile_options(
    ${PROJECT}
    PRIVATE
    $<$<CONFIG:Release>:-mfma>
    $<$<CONFIG:Debug>:-mfma>
    )
endif()

enable_testing()
add_subdirectory(utest)
