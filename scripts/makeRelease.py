#!/usr/bin/python3



"""

This script creates a tar file containing a release package containing 
shared libraries for all supported platforms.

Is should be invoked with a single argument, the name of the release 
package to create (e. g. 0.1.0). 

This script must be invoked from a directory containing the ExpressionMatrix2 tree.
The following goes into the tar file:
ExpressionMatrix2/doc
ExpressionMatrix2/tests
ExpressionMatrix2/platformName/ExpressionMatrix2.so (for all supported platforms).
The shared libraries are assumed to be up to date. This script
does not attempt to rebuild them. Remember to strip each of them
before creating the release package!

"""



# Import what we need.
import os
import sys



# Function to check if a path exists and is a regular file or a directory.
def regularFileExists(name):
    return os.path.exists(name) and os.path.isfile(name)
def directoryExists(name):
    return os.path.exists(name) and os.path.isdir(name)



# Get the releaseName as the first and only argument.
if len(sys.argv) != 2:
    raise Exception('Invoke with one parameter, the name of the release (e. g. 0.0.1).')
releaseName = sys.argv[1]



# Make sure we are in the right directory.
if (
    not directoryExists('ExpressionMatrix2') or 
    not directoryExists('ExpressionMatrix2/src') or 
    not directoryExists('ExpressionMatrix2/doc') or 
    not directoryExists('ExpressionMatrix2/tests')):
    raise Exception('This command should be invoked from the directory containing the ExpressionMatrix2 directory.')



# Create the list of files and directories to add to the tar file.
files = ['ExpressionMatrix2/doc', 'ExpressionMatrix2/tests']
platformNames = [
    'Release-ubuntu16-python3',
    'Release-centos7-python2',
]
for platformName in platformNames:
    platformDirectory = 'ExpressionMatrix2/' + platformName

    # Check that the directory for this platform exists.
    if not directoryExists(platformDirectory):
        raise Exception('The platform directory ' + platformDirectory + ' does not exist.')

    # Make sure the shared library exists. 
    # We assume that is is up to date and not attempt to rebuild it.
    sharedLibraryName = platformDirectory + '/ExpressionMatrix2.so'   
    if not regularFileExists(sharedLibraryName):
        raise Exception(sharedLibraryName + ' not found')

    # Add it to the list of files to be included in the tar file.
    files.append(sharedLibraryName)


# Create the tar file.
tarFileName = 'ExpressionMatrix2-' + releaseName + '.tar'
os.system('tar --create --bzip2 --file ' + tarFileName + ' ' + ' '.join(files))






