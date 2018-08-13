cmake_minimum_required(VERSION 3.2)

set(PROJECT plbpd)
project(${PROJECT} CXX)

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, defaulting to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

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
option(AVX2 "Build plbpd with AVX2 intrinsics" ON)
if(AVX2)
  message(STATUS "Adding -DAVX2...")
  add_definitions(-DAVX2)
  target_compile_options(
    ${PROJECT}
    PRIVATE
    $<$<CONFIG:Release>:-mavx2 -march=core-avx2>
    $<$<CONFIG:Debug>:-mavx2>
    $<$<CONFIG:RelWithDebInfo>:-mavx2>
    $<$<CONFIG:MinSizeRel>:-mavx2>
    )
endif()

enable_testing()
add_subdirectory(utest)
