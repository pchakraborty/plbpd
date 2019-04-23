cmake_minimum_required(VERSION 3.2)
set (CMAKE_CXX_STANDARD 14)

add_executable(TestVectorField TestVectorField.cpp)
add_test(NAME EqlbDist COMMAND TestVectorField)
