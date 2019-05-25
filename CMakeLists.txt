cmake_minimum_required(VERSION 3.10)

project(cm_heat_equation
        DESCRIPTION "Heat equation solver tool."
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

set(sources main.cpp)
set(cli CLI/CLI11.hpp)

add_executable(cm_heat_equation ${sources} ${cli})

add_subdirectory(heateq)
target_link_libraries(cm_heat_equation heat_equation)

target_compile_options(cm_heat_equation PRIVATE
        $<$<CONFIG:DEBUG>: -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC>
        $<$<CONFIG:RELEASE>: -DNDEBUG -O2 -funroll-loops -ftree-vectorize>)
target_compile_options(cm_heat_equation PRIVATE $<$<CXX_COMPILER_ID:GNU>:
        $<$<CONFIG:DEBUG>: -ggdb -Wformat=2 -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align>>)
target_compile_options(cm_heat_equation PRIVATE $<$<CXX_COMPILER_ID:Clang>:
        $<$<CONFIG:DEBUG>:-Weverything -Wno-c++98-compat-pedantic -Wno-missing-prototypes>>)
