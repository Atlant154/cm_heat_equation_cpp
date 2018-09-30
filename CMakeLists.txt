cmake_minimum_required(VERSION 3.10)
project(cm_heat_equation)

set(CMAKE_CXX_STANDARD 14)

set(sources src/main.cpp src/tridiagonal_matrix.cpp)
set(headers include/tridiagonal_matrix.h)

add_executable(vychi ${sources} ${headers})

if(GNU)
	set(CMAKE_CXX_FLAGS_DEBUG  "-O2 -Wall -Wextra -Wshadow -pedantic")
	set(CMAKE_CXX_FLAGS_RELEASE  "-O3")
endif()