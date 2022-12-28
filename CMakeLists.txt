cmake_minimum_required(VERSION 3.23)
project(lab_3)

set(CMAKE_CXX_STANDARD 20)

add_executable(lab_3 main.cpp graph.cpp graph.h)

add_compile_options(-O2)
add_compile_options(-march=native)
