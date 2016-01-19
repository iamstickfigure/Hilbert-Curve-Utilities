cmake_minimum_required(VERSION 3.1)
project(Hilbert-Curve-Utilities)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
link_directories("/usr/local/include")

set(SOURCE_FILES cpp/hilbert_util.cpp)
add_executable(untitled ${SOURCE_FILES})