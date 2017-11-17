#!/usr/bin/python3



"""

This script creates a tar file containing a release package containing 
shared libraries for all supported platforms.

Is should be invoked with a single argument, the name of the release 
package to create (e. g. 0.1.0). 

This script must be invoked from the top level directory of the 
ExpressionMatrix2 tree (that is, from the directory containing the src, doc, 
and tests directories).
The following goes into the tar file:
doc
tests
platformName/ExpressionMatrix2.so (for all supported platforms).
The shared libraries are assumed to be up to date. This script
does not attempt to rebuild them. 



PROCEDURE TO CREATE A NEW RELEASE

- Commit and push all pending code.
- For each of the supported platforms:
  * Get the latest code.
  * Build ExpressionMatrix2.so.
  * Strip ExpressionMatrix2.so.
  * Check that all the tests work.
  * Copy ExpressionMatrix2.so to the ExpressionMatrix2 tree where makeRelease.py will run.
- Create updated Python documentation and push it to GitHub.
- Run makeRelease.py passing the release number (e. g. 0.1.0) as the argument.
  This creates the tar file for the release.
- Verify that the size of the tar file is reasonable.
- Verify that the tar file extracts correctly to a directory ExpressionMatrix2-x.y.z,
  where x.y.z is the release name.
- In GitHub:
  * Create the new release. 
  * Upload the tar file.
  * Write some release notes.
- Remove the local copy of the tar file.

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
    not directoryExists('src') or 
    not directoryExists('doc') or 
    not directoryExists('tests')):
    raise Exception('This command should be invoked from the top level of the ExpressionMatrix2 tree.')



# Create the list of files and directories to add to the tar file.
files = ['doc', 'tests']
platformNames = [
    'Release-ubuntu16-python3',
    'Release-centos7-python2',
]
for platformName in platformNames:

    # Check that the directory for this platform exists.
    if not directoryExists(platformName):
        raise Exception('The platform directory ' + platformName + ' does not exist.')

    # Make sure the shared library exists. 
    # We assume that is is up to date and not attempt to rebuild it.
    sharedLibraryName = platformName + '/ExpressionMatrix2.so'   
    if not regularFileExists(sharedLibraryName):
        raise Exception(sharedLibraryName + ' not found')

    # Add it to the list of files to be included in the tar file.
    files.append(sharedLibraryName)


# Create the tar file.
tarFileName = 'ExpressionMatrix2-' + releaseName + '.tar'
transformOption = "--transform 's,^,ExpressionMatrix2-" + releaseName + "/,'"
os.system(' '.join(['tar --create --bzip2 --file', tarFileName, transformOption] + files))






