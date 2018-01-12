#!/usr/bin/python3

"""

This script runs the toy tests on all the platforms
build by the makeAll script appropriate for the
system we are running on.

"""

# Import what we need.
import glob
import os
import shutil
import sys


# Get the path of the code tree to test.
if len(sys.argv) != 2:
    raise Exception('Invoke with one argument, the path to the code tree to be tested.')
path = sys.argv[1]

# Check that it exists.
if not os.path.exists(path):
    raise Exception(path + ' does not exist.')    
    
# cd to the build directory.
os.chdir(path)
os.chdir('ExpressionMatrix2/build')



# Loop over all platforms we find.
for platform in glob.glob('*'):
    print('Testing ' + platform)
    
    # Find the name of the python command for this platform
    if platform.find('python3') >= 0:
        pythonCommand = 'python3'
    else:
        pythonCommand = 'python'
    print(pythonCommand)
    
    # Go to the build directory for this platform.
    os.chdir(platform)
    
    # Point to the ExpressionMatrix2.so here.
    os.environ['PYTHONPATH'] = os.getcwd()
    
    # Copy the toy tests here.
    os.system('cp -rp ../../tests/ToyTest* .') 
    
    # Loop over the toy tests.
    for test in glob.glob('ToyTest*'):
        print('Testing %s on %s.' % (test, platform))
        os.chdir(test)
        print('Test directory is ' + os.getcwd())
        if os.path.exists('data'):
            shutil.rmtree('data')
        returnCode = os.system(pythonCommand + ' run.py')
        if returnCode:
            raise Exception('Test %s failed on %s.' % (test, platform))
        os.chdir('..')

    # Go back to the build directory.
    os.chdir('..')
    




