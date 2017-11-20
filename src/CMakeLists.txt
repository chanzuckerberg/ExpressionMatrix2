# Cmake file for ExpressionMatrix2/src.

# This creates shared library ExpressionMatrix2.so,
# which can be imported as a Python package.

# To build:
# - Create a new empty directory to contain the build and cd to it.
#   This directory can be anywhere - it does not have to
#   be in the ExpressionMatrix2 source tree.
# - Run command "cmake .../ExpressionMatrix2/src",
#   making sure to enter the local path for the ExpressionMatrix2/src
#   directory. If necessary, use cmake option -D to override
#   default values for one or more of the following 
#   configuration variables, listed here with their
#   default values suitable for ubuntu16 and Python 3:
#   BUILD_WITH_HDF5=ON
#   PYTHON_INCLUDE_PATH=/usr/include/python3.5m
#   PYBIND11_INCLUDE_PATH=/usr/local/include/python3.5
#   HDF5_INCLUDE_PATH=/usr/include/hdf5/serial
#   HDF5_LIBRARIES=hdf5_cpp hdf5_serial
# - Run command "make". If your machine has more than
#   one processor and sufficient memory you can speed up to build
#   using the make "-j" option - see the man page for the make
#   command for details.

#   Here are examples of -D options that work for various platforms:
#   Ubuntu 16, Python 3: 
#       No option needed.
#   Ubuntu 16, Python 2: 
#       -DPYTHON_INCLUDE_PATH=/usr/include/python2.7 
#       -DPYBIND11_INCLUDE_PATH=/usr/local/include/python2.7
  
# This file is only appropriate for a release build on Linux.
# The resulting file can be sigificantly reduced in size
# using the strip command. This cmake file
# does not do that, to create a release build
# that is debuggable.

cmake_minimum_required(VERSION 2.6)
project(ExpressionMatrix2)


# Create a Makefile that writes commands before executing them.
SET(CMAKE_VERBOSE_MAKEFILE ON)

# The source files.
# Done like this, cmake needs to rerun when new files are added.
file(GLOB SOURCES *.cpp)

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

add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)



# Define our shared library.
add_library(ExpressionMatrix2 SHARED ${SOURCES})

# To make sure the library is named ExpressionMatrix2.so,
# get rid of the "lib" prefix.
set_target_properties(ExpressionMatrix2 PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(ExpressionMatrix2 PROPERTIES  DEFINE_SYMBOL "")

# Boost libraries.
# target_link_libraries(ExpressionMatrix2 boost_filesystem)


# Definition of configuration variables begins here.
message(STATUS 
    "\nConfiguration variables: use cmake option -D to override:")



# Option to turn off hdf5 functionality (to eliminate
# dependency on hdf5 include files and libraries).
option(BUILD_WITH_HDF5 "Include HDF5 functionality in the build." ON)
message(STATUS "BUILD_WITH_HDF5=" ${BUILD_WITH_HDF5})
if(NOT BUILD_WITH_HDF5)
    add_definitions(-DCZI_EXPRESSION_MATRIX2_SKIP_HDF5)
endif(NOT BUILD_WITH_HDF5)


# Include directory for Python.
# This determines the Python version that the
# library will work with.
if(NOT PYTHON_INCLUDE_PATH)
    # Default value for ubuntu16 and Python 3.
    set(PYTHON_INCLUDE_PATH /usr/include/python3.5m)
endif(NOT PYTHON_INCLUDE_PATH)
include_directories(${PYTHON_INCLUDE_PATH})
message(STATUS "PYTHON_INCLUDE_PATH=" ${PYTHON_INCLUDE_PATH})



# Include directory for pybind11.
if(NOT PYBIND11_INCLUDE_PATH)
    # Default value for ubuntu16 and Python 3.
    set(PYBIND11_INCLUDE_PATH /usr/local/include/python3.5)
endif(NOT PYBIND11_INCLUDE_PATH)
include_directories(${PYBIND11_INCLUDE_PATH})
message(STATUS "PYBIND11_INCLUDE_PATH=" ${PYBIND11_INCLUDE_PATH})



# Include directory for hdf5.
if(NOT HDF5_INCLUDE_PATH)
    # Default value for ubuntu16.
    set(HDF5_INCLUDE_PATH /usr/include/hdf5/serial)
endif(NOT HDF5_INCLUDE_PATH)
if(BUILD_WITH_HDF5)
    include_directories(${HDF5_INCLUDE_PATH})
endif(BUILD_WITH_HDF5)
message(STATUS "HDF5_INCLUDE_PATH=" ${HDF5_INCLUDE_PATH})



# Hdf5 libraries.
if(NOT HDF5_LIBRARIES)
    # Default value for ubuntu16.
    set(HDF5_LIBRARIES hdf5_cpp hdf5_serial)
else(NOT HDF5_LIBRARIES)
    separate_arguments(HDF5_LIBRARIES)
endif(NOT HDF5_LIBRARIES)
string(REPLACE ";" " " SPACE_SEPARATED_HDF5_LIBRARIES "${HDF5_LIBRARIES}")
message(STATUS HDF5_LIBRARIES=${SPACE_SEPARATED_HDF5_LIBRARIES})
if(BUILD_WITH_HDF5)
    foreach(HDF5_LIBRARY IN LISTS HDF5_LIBRARIES)
        target_link_libraries(ExpressionMatrix2 ${HDF5_LIBRARY})
    endforeach(HDF5_LIBRARY)
endif(BUILD_WITH_HDF5)



# Message that concludes the list of configuration variables.
message(STATUS "\n")


