cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(watcor)

add_executable(watcor main.cpp readall.cpp gro.cpp model.cpp)

install(TARGETS watcor RUNTIME DESTINATION bin)
