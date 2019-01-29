cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
project(heat_equation_tests)

include(GoogleTest)
add_subdirectory(google_test EXCLUDE_FROM_ALL)

set(gtest_force_shared_crt OFF CACHE BOOL "" FORCE)

macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main heat_equation)
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

package_add_test(AccuracyTest accuracy.cpp)