cmake_minimum_required(VERSION 3.10)
project(cm_heat_equation)

set(CMAKE_CXX_STANDARD 17)

set(sources main.cpp src/heat_equation.cpp)
set(headers include/heat_equation.h)

add_executable(cm_heat_equation ${sources} ${headers})

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra -pedantic -g -Wshadow -Wformat=2 -Wconversion \
        -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC \
        -D_FORTIFY_SOURCE=2 -fsanitize=undefined")

    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -fno-tree-vectorize -fprofile-use -ftree-loop-vectorize")

endif ()