cmake_minimum_required(VERSION 3.5)
project(new_rkf45_with_lens)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp BaseLens.cpp BaseLens.h Lens.h Lens.cpp Coords.h lightRay.cpp lightRay.h)
add_executable(new_rkf45_with_lens ${SOURCE_FILES})