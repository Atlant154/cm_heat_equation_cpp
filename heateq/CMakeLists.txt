cmake_minimum_required(VERSION 3.10)

project(heat_equation
        VERSION 2.0
        DESCRIPTION "Numeric solution of heat equation ODE library."
        LANGUAGES CXX
)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui:
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

set( headers include/HeatEquation.hpp serialization/heq.pb.h )
set( sources src/HeatEquation.cpp serialization/heq.pb.cc )

add_library(heat_equation STATIC ${headers} ${sources})
target_include_directories(heat_equation PUBLIC include)
target_include_directories(heat_equation PRIVATE serialization)

find_package(Protobuf REQUIRED)

target_link_libraries(heat_equation tbb ${Protobuf_LIBRARIES})

set_target_properties(heat_equation PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS OFF
        INTERPROCEDURAL_OPTIMIZATION True
)


find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CUDA_COMPILER_LAUNCHER "${CCACHE_PROGRAM}") # CMake 3.9+
endif()


target_compile_options(heat_equation PRIVATE
        $<$<CONFIG:DEBUG>: -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC>
        $<$<CONFIG:Release>: -DNDEBUG -O2 -funroll-loops -ftree-vectorize>
)

target_compile_options(heat_equation PRIVATE $<$<CXX_COMPILER_ID:GNU>:
        $<$<CONFIG:DEBUG>: -ggdb -Wformat=2 -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align>>
)

target_compile_options(heat_equation PRIVATE $<$<CXX_COMPILER_ID:Clang>:
        $<$<CONFIG:DEBUG>:-Weverything -Wno-c++98-compat-pedantic -Wno-missing-prototypes>>
)

add_subdirectory(benchmark)
add_subdirectory(tests)