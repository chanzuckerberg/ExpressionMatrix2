cmake_minimum_required(VERSION 2.6)
project(ExpressionMatrix2)

# Cmake file for ExpressionMatrix2/src.
# This creates shared library ExpressionMatrix2.so,
# which can be imported as a Python package.

# To build:
# - Create a new empty directory to contain the build and cd to it.
#   This directory can be anywhere - it does not have to
#   be in the ExpressionMatrix2 source tree.
# - Run command "cmake .../ExpressionMatrix2/src",
#   making sure to enter the local path for the ExpressionMatrix2/src
#   directory.
# - Run command "make". If your machine has more than
#   one processor and sufficient memory you can speed up to build
#   using the make "-j" option - see the man page for the make
#   command for details.
  
# This file is only appropriate for a release build on Linux.
# The resulting file can be sigificantly reduced in size
# using the strip command. This cmake file
# does not do that, to create a release build
# that is debuggable.



# The source files.
# Done like this, cmake needs to rerun when new files are added.
file(GLOB SOURCES *.cpp *.hpp)

# Compile flag to specify the C++ standard we are using.
# This corresponds to C+11.
# Newer versions of cmake have special provisions
# for this, but we don't want to rely on that.
add_definitions(-std=c++0x)

# Generate debug information in gdb format.
# This way we create a debuggable release build.
# If this is not needed, uou can use the strip command
# to drastically reduce the size of the library.
add_definitions(-ggdb3)

# Compile flags to control optimization.
# These are appropriate for a release build.
# The -msse4.2 flag is necessary to enable use
# of the popcount instruction, which is essential
# for performance of the LSH code.
add_definitions(-O3 -msse4.2)

# Options to control compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result)

# Include directory for Python.
include_directories(/usr/include/python3.5m)

# Include directory for pybind11.
include_directories(/usr/local/include/python3.5)

# Include directory for hdf5.
include_directories(/usr/include/hdf5/serial)

# Create a Makefile that writes commands before executing them.
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Define our shared library.
add_library(ExpressionMatrix2 SHARED ${SOURCES})

# To make sure the library is named ExpressionMatrix2.so,
# get rid of the "lib" prefix.
set_target_properties(ExpressionMatrix2 PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(ExpressionMatrix2 PROPERTIES  DEFINE_SYMBOL "")

# Boost libraries.
target_link_libraries(ExpressionMatrix2 boost_filesystem)

# HDF5 libraries.
target_link_libraries(ExpressionMatrix2 hdf5_cpp hdf5_serial)

