cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
project(heat_equation DESCRIPTION "Numeric solution of heat equation ODE library.")

add_library(heat_equation SHARED src/HeatEquation.cpp)
target_include_directories(heat_equation PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include/>
        PRIVATE ${CMAKE_CURRENT_LIST_DIR}/src)

add_subdirectory(benchmark)

target_compile_options(cm_heat_equation PRIVATE
        $<$<CONFIG:DEBUG>: -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC>
        $<$<CONFIG:RELEASE>: -DNDEBUG -march=native -O2 -funroll-loops -ftree-vectorize>)
target_compile_options(cm_heat_equation PRIVATE $<$<CXX_COMPILER_ID:GNU>:
        $<$<CONFIG:DEBUG>: -ggdb -Wformat=2 -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align>>)
target_compile_options(cm_heat_equation PRIVATE $<$<CXX_COMPILER_ID:Clang>:
        $<$<CONFIG:DEBUG>:-Weverything -Wno-c++98-compat-pedantic -Wno-missing-prototypes>>)
target_link_libraries(heat_equation stdc++fs)
