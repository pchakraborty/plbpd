cmake_minimum_required(VERSION 3.2)

add_executable(TestEqlbDist ../LBModel.cpp TestEqlbDist.cpp)
add_test(NAME EqlbDist COMMAND TestEqlbDist)
