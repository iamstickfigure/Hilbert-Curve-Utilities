cmake_minimum_required(VERSION 3.1)
project(Hilbert-Curve-Utilities)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories("/usr/local/include")
#add_library(libgraph /usr/local/include/graphics.h)

set(SOURCE_FILES cpp/hilbert_util.cpp)
add_executable(hilbert_util cpp/hilbert_util.cpp)
#target_link_libraries(hilbert_util libgraph)