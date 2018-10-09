cmake_minimum_required(VERSION 3.10)
project(cm_heat_equation)

set(CMAKE_CXX_STANDARD 14)

set(sources src/main.cpp src/tridiagonal_matrix.cpp)
set(headers include/tridiagonal_matrix.h)

add_executable(cm_heat_equation ${sources} ${headers})

if(GNU)
    set(CMAKE_CXX_FLAGS_DEBUG  "-O2 -Wall -Wextra -pedantic -Wshadow -Wformat=2 -Wconversion \
                                -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align 
                                -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC \
                                -D_FORTIFY_SOURCE=2 -fstack-protector")
    set(CMAKE_CXX_FLAGS_RELEASE  "-O3")
endif()