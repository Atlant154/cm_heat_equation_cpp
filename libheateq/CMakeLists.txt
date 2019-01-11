cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
project(heat_equation DESCRIPTION "Numeric solution of heat equation ODE library.")

add_library(heat_equation SHARED src/heat_equation.cpp)
set_target_properties(heat_equation PROPERTIES PUBLIC_HEADER include/heat_equation.h)

add_subdirectory(benchmark)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra -pedantic -g -Wshadow -Wformat=2 -Wconversion \
        -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC \
        -D_FORTIFY_SOURCE=2 -fsanitize=undefined")

    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -funroll-loops -ftree-vectorize")

endif ()