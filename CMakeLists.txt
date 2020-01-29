cmake_minimum_required(VERSION 3.15)
project(radLength)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR})

add_executable(radLength main.cpp cxxopts.hpp)