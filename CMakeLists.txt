cmake_minimum_required(VERSION 3.10)
project(cm_heat_equation)

set(CMAKE_CXX_STANDARD 17)

set(sources main.cpp)
add_executable(cm_heat_equation ${sources})

add_subdirectory(heateq)
target_link_libraries(cm_heat_equation heat_equation)
enable_testing()

option(CLI11_TESTING "Build the tests and add them" OFF)
add_subdirectory(CLI11)
include_directories(CLI11/include)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra -pedantic -g -Wshadow -Wformat=2 -Wconversion \
        -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC \
        -D_FORTIFY_SOURCE=2 -fsanitize=undefined")

    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -funroll-loops -ftree-vectorize")

endif ()
