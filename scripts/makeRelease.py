#!/usr/bin/python3

import os

"""
This scripts creates a tar file containing a release package for the current, existing build.
It should be invoked fom the top level directory of the ExpressionMatrix2 tree
(that is, the directory containing src, doc, etc.). It will create a file named
Release.tar in the same directory.
"""

# Function to check if a directory exists.
def directoryExists(name):
    return os.path.exists(name) and os.path.isdir(name)

# Make sure we are in the right directory.
message = 'The src directory is missing. This command should be invoked from the top level directory of the ExpressionMatrix2 directory.'
currentDirectory = os.getcwd()
if not directoryExists('src') or not directoryExists('doc') or not directoryExists('Release') or not directoryExists('tests'):
    raise Exception(message)

# Make sure the shared library exists.    
if not os.path.exists('Release/ExpressionMatrix2.so'):
    raise Exception('Release/ExpressionMatrix2.so not found')

# Create the tar file with the doc and tests directories.
os.system('tar --create --file Release.tar doc tests')

# Strip the shared library, then add it to the archive at the top level.
# Stripping it reduces its size by a large factor.    
os.chdir('Release')  
os.system('strip ExpressionMatrix2.so')
os.system('tar --append --file ../Release.tar ExpressionMatrix2.so')


