# Cmake file for directory ExpressionMatrix2/src-static.

# This builds executable ExpressionMatrix2,
# which can be useful in situations where the Python API
# is not needed and the user just wants to use the http server.

# This executable is statically linked and so it 
# does not depend on any shared library
# and is highly portable to most modern Linux distributions.

# On 3/16/2018, the executable built on Ubuntu 16.04
# ran without problems on CentOS 7, ArchLinux, 
# Fedora 27, and Ubuntu 18.04 Beta.

# For complete functionality, the graphviz package is required
# as it is used to compute graph layouts.
# The graphviz package is easily available for most distributions.

# This executable is called without arguments.
# It creates or accesses the ExpressionMatrix object
# stored in directory data, then starts the http server on port 17100.

cmake_minimum_required(VERSION 2.6)
project(ExpressionMatrix2-static)

file(GLOB SOURCES main.cpp ../src/*.cpp)
include_directories(../src)

add_executable(ExpressionMatrix2 ${SOURCES})
set_target_properties(ExpressionMatrix2 PROPERTIES LINK_FLAGS "-static" )

add_definitions(-std=c++0x)
add_definitions(-ggdb3)
add_definitions(-O3 -msse4.2)
add_definitions(-Wall -Wconversion -Wno-unused-result)

add_definitions(-DCZI_EXPRESSION_MATRIX2_SKIP_HDF5)
target_link_libraries(ExpressionMatrix2 pthread)

SET(CMAKE_VERBOSE_MAKEFILE ON)



