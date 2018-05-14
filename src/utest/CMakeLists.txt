cmake_minimum_required(VERSION 3.2)

add_executable(TestEqlbDist ../LBModel.cpp TestEqlbDist.cpp)
add_test(NAME EqlbDist COMMAND TestEqlbDist)

add_executable(TestLBModel ../LBModel.cpp TestLBModel.cpp)
add_test(NAME LBModel COMMAND TestLBModel)

add_executable(TestDomain ../LBModel.cpp ../Domain.cpp TestDomain.cpp)
add_test(NAME Domain COMMAND TestDomain)

add_executable(TestArrayND TestArrayND.cpp)
add_test(NAME ArrayND COMMAND TestArrayND)
