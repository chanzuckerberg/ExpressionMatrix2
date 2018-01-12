#!/usr/bin/python3

"""

This script can be used to create a tar file to be uploaded
to GitHub for a release.

This script should run on the Ubuntu 16 build machine
passing two arguments:
- The name of the release to be created (for example, 0.4.0).
- The ssh argument (userName@machineName) for the remote CentOS 7 machine.
  It is assumed that public/private keys are properly set up
  so it is possible to issue ssh commands without entering a password.
  
This clones the repository and does a build and test on both 
the local Ubuntu 16 machine and the remote CentOS 7 machine.
On each of the machines, all results go to directory
/tmp/ExpressionMatrix2-ReleaseName, which should not already exist.
Build output on the remote CentOS 7 machine are then copied to
the local Ubuntu 16 machine, and a tar file is generated
as /tmp/ExpressionMatrix2-ReleaseName/ExpressionMatrix2/ExpressionMatrix2-ReleaseName.tar.

The files created are not removed. They should be removed manually if desired.



PROCEDURE TO CREATE A NEW RELEASE

- Make sure the documentation for the Python API is up to date. 
- Run this script on the local Ubuntu 16 machine to create the tar file.
- Verify that the size of the tar file is reasonable.
- Verify that the tar file extracts correctly to a directory ExpressionMatrix2-ReleaseName.
- In GitHub:
  * Create the new release. 
  * Upload the tar file.
  * Write some release notes.
  * Publish the release.
- If desired, remove the files created in /tmp/ExpressionMatrix2-ReleaseName
  on both the local Ubuntu 16 machine and the remote CentOS 7 machine.

"""

# Import what we need.
import glob
import os
import sys

# Get the arguments.
if len(sys.argv) != 3:
    raise Exception(
    'Invoke with two arguments:\n' 
    '- The name of the release to be created.\n'
    '- The ssh argument (userName@machineName) for the remote CentOS 7 machine.\n')
releaseName = sys.argv[1]
centos7Machine = sys.argv[2]
workDirectory = '/tmp/ExpressionMatrix2-' + releaseName;

# Locate the scripts we need.
makeAllScript = os.path.dirname(__file__) + '/makeAll.py'
runToyTests = os.path.dirname(__file__) + '/runToyTests.py'

# Run makeAll.py on the ubuntu16 7 machine.
os.system('%s ubuntu16 %s' % (makeAllScript, workDirectory))

# Run runToyTests.py on the ubuntu16 7 machine.
os.system('%s %s' % (runToyTests, workDirectory))

# Run makeAll.py on the CentOS 7 machine.
returnCode = os.system(' ssh %s python < %s - centos7 %s' % (centos7Machine, makeAllScript, workDirectory))
if returnCode:
    raise Exception('CentOS 7 build failed.')

# Run runToyTests.py on the CentOS 7 machine.
returnCode = os.system(' ssh %s python < %s - %s' % (centos7Machine, runToyTests, workDirectory))
if returnCode:
    raise Exception('CentOS 7 tests failed.')
    
# Copy the build directory of the CentOS 7 machine to the local (ubuntu16) machine.
os.system('scp -rp %s:%s/ExpressionMatrix2/build %s/ExpressionMatrix2/build-centos7' % 
    (centos7Machine, workDirectory, workDirectory)) 
os.system('mv %s/ExpressionMatrix2/build-centos7/* %s/ExpressionMatrix2/build' % (workDirectory, workDirectory))
os.system('rmdir %s/ExpressionMatrix2/build-centos7' % workDirectory)



# At this point all the files we need to create the tar file are available locally.
os.chdir('%s/ExpressionMatrix2' % workDirectory)

# Create the list of files and directories to add to the tar file.
files = ['doc', 'tests'] + glob.glob('build/*/ExpressionMatrix2.so')

# Create the tar file.
tarFileName = 'ExpressionMatrix2-' + releaseName + '.tar'
transformOption = "--transform 's,^,ExpressionMatrix2-" + releaseName + "/,'"
os.system(' '.join(['tar --create --bzip2 --file', tarFileName, transformOption] + files))
print('Created tar file %s/ExpressionMatrix2/%s' % (workDirectory, tarFileName))


