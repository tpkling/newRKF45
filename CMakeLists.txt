cmake_minimum_required(VERSION 3.5)
project(new_rkf45_with_lens)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp Lens.h Lens.cpp BaseLens.cpp BaseLens.h randomGen.cpp randomGen.h LensGroup.cpp LensGroup.h lightRay.cpp lightRay.h SourceFuncs.cpp SourceFuncs.h)
add_executable(new_rkf45_with_lens ${SOURCE_FILES})