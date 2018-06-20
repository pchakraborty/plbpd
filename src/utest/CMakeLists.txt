cmake_minimum_required(VERSION 3.2)

add_executable(TestEqlbDist ../LBModel.cpp TestEqlbDist.cpp)
add_test(NAME EqlbDist COMMAND TestEqlbDist)

add_executable(TestLBModel ../LBModel.cpp TestLBModel.cpp)
add_test(NAME LBModel COMMAND TestLBModel)

add_executable(TestDomain ../LBModel.cpp ../Domain.cpp TestDomain.cpp)
add_test(NAME Domain COMMAND TestDomain)

add_executable(TestField TestField.cpp)
add_test(NAME Field COMMAND TestField)

include(../cmake/hdf5.cmake)
include_directories(${HDF5_INCLUDE_DIR})
add_executable(TestLattice ../Lattice.cpp TestLattice.cpp)
target_link_libraries(TestLattice ${HDF5_C_STATIC_LIBRARY})
add_test(NAME Lattice COMMAND TestLattice)

add_test(NAME plbpd COMMAND plbpd)
  # NAME plbpd
  # COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/plbpd.cmake
  # )
