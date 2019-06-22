cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
project(heat_equation_tests)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} ${GTEST_LIBRARIES} pthread heat_equation)
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

package_add_test(AccuracyTest accuracy.cpp)