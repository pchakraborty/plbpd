cmake_minimum_required(VERSION 3.2)

add_executable(TestEqlbDist ../LBModel.cpp TestEqlbDist.cpp)
add_test(NAME EqlbDist COMMAND TestEqlbDist)

add_executable(TestDomain ../LBModel.cpp ../Domain.cpp TestDomain.cpp)
add_test(NAME Domain COMMAND TestDomain)
