#!/usr/bin/python3

"""

This script does a clone from GitHub to create 
a brand new code tree at the specified path
and builds all of the supported configurations for the
system we are running on.

Run with two arguments:
- The system we are running on.
  Currently, this can only be ubuntu16 or centos7.
- The path to the new code tree to be created.
It must not exist.

This scripts assumes that all packages required by the build
are installed.

"""

# Import what we need.
import os
import sys

# Check the number of arguments.
if len(sys.argv) != 3:
    raise Exception(
    'Invoke with two arguments:\n' 
    '- The system to build for (ubuntu16 or centos7).\n'
    '- The path to the code tree to be created.\n')


# The platforms their options.
ubuntu16Platforms = {
    'Release-ubuntu16-python3' :
    [
    ],
    'Release-ubuntu16-python2' :
    [
    '-DPYTHON_INCLUDE_PATH=/usr/include/python2.7', 
    '-DPYBIND11_INCLUDE_PATH=/usr/local/include/python2.7',
    ],
    'Release-ubuntu16-nohdf5-python3' :
    [
    '-DBUILD_WITH_HDF5=OFF',
    ],
    'Release-ubuntu16-nohdf5-python2' :
    [
    '-DPYTHON_INCLUDE_PATH=/usr/include/python2.7', 
    '-DPYBIND11_INCLUDE_PATH=/usr/local/include/python2.7',
    '-DBUILD_WITH_HDF5=OFF',
    ],
    }
centos7Platforms = {
    'Release-centos7-python2' :
    [
    '-DPYTHON_INCLUDE_PATH=/usr/include/python2.7', 
    '-DPYBIND11_INCLUDE_PATH=/usr/lib/python2.7/site-packages',
    '-DHDF5_INCLUDE_PATH=/usr/include',
    '-DHDF5_LIBRARIES="hdf5_cpp hdf5"',
    ],
    'Release-centos7-nohdf5-python2' :
    [
    '-DPYTHON_INCLUDE_PATH=/usr/include/python2.7', 
    '-DPYBIND11_INCLUDE_PATH=/usr/lib/python2.7/site-packages',
    '-DBUILD_WITH_HDF5=OFF',
    ],
    }

allPlatforms = {
    'ubuntu16': ubuntu16Platforms,
    'centos7':  centos7Platforms,
    }

# Get the platforms we are going to build
systemName = sys.argv[1]
if not systemName in allPlatforms.keys():
    raise Exception('The first argument must be one of ' + ' '.join(allPlatforms.keys()))
platforms = allPlatforms[systemName]
    

# Get the path where the new code tree should be created.
path = sys.argv[2]

# Check that it does not already exists.
if os.path.exists(path):
    raise Exception(path + ' already exists.')    
    
# Create it and cd to it.
os.mkdir(path)
os.chdir(path)

# Clone the GitHub repository.
os.system('git clone https://github.com/chanzuckerberg/ExpressionMatrix2')
os.chdir('ExpressionMatrix2')

# Create the build directory and cd to it.
os.mkdir('build')
os.chdir('build')

# Build all platforms
for platform in platforms:
    print('Building ' + platform)
    options = platforms[platform]
    os.mkdir(platform)
    os.chdir(platform)
    returnCode = os.system('cmake ../../src ' + ' '.join(options))
    if returnCode:
        raise Exception('cmake error for ' + platform)
    returnCode = os.system('make -j all')
    if returnCode:
        raise Exception('Build error for ' + platform)
    os.system('strip ExpressionMatrix2.so')
    os.chdir('..')





