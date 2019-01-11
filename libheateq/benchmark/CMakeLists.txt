cmake_minimum_required(VERSION 3.10)

project(heat_equation_benchmark)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
option(BENCHMARK_DOWNLOAD_DEPENDENCIES "Download gBenchmark gTest deps" ON)

add_subdirectory(google_benchmark)
include_directories(google_benchmark/include)
add_executable(heat_equation_benchmark benchmark.cpp)

target_link_libraries(heat_equation_benchmark benchmark)
target_link_libraries(heat_equation_benchmark heat_equation)
