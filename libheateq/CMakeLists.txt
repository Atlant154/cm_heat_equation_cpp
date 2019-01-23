cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
project(heat_equation DESCRIPTION "Numeric solution of heat equation ODE library.")

add_library(heat_equation SHARED src/heat_equation.cpp)
target_include_directories(heat_equation PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include/>
        PRIVATE ${CMAKE_CURRENT_LIST_DIR}/src)

add_subdirectory(benchmark)

add_subdirectory(json)
target_link_libraries(heat_equation nlohmann_json::nlohmann_json)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra -pedantic -g -Wshadow -Wformat=2 -Wconversion \
        -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC \
        -D_FORTIFY_SOURCE=2 -fsanitize=undefined")

    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -funroll-loops -ftree-vectorize")

    target_link_libraries(heat_equation stdc++fs)

endif ()

option(PACKAGE_TESTS "Build the tests" ON)
if (PACKAGE_TESTS)

    if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    endif ()

    enable_testing()
    add_subdirectory(tests)
endif ()
