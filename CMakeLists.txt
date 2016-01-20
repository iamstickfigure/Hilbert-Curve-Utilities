cmake_minimum_required(VERSION 3.1)
project(Hilbert-Curve-Utilities)

set(GCC_CIMG_COMPILE_FLAGS "-O2 -L/usr/X11R6/lib -lm -lpthread -lX11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${GCC_CIMG_COMPILE_FLAGS}")
#include_directories("/usr/local/include")
include_directories("/usr/local/include /usr/include")

set(SOURCE_FILES cpp/hilbert_util.cpp)
add_executable(hilbert_util cpp/hilbert_util.cpp)