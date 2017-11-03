# ExpressionMatrix2

This repository contains software for analysis, visualization, and clustering of gene expression data from single-cell RNA sequencing developed at [Chan-Zuckerberg Initiative](https://chanzuckerberg.com).
It scales favorably to large numbers of cells thank to its use of Locality-Sensitive Hashing (LSH), and was successfully used, without downsampling, on a data set with over one million cells.

Documentation for the latest version of this software is available online through [GitHub Pages](https://chanzuckerberg.github.io/ExpressionMatrix2/doc/index.html), or you can use the directions below to obtain documentation for any previous release.  

This code is at a prototype, pre-alpha stage. It is hoped that it can be useful in its current form, but it is likely that it contains bugs or other errors that impact its ability to give meaningful results. In addition, a validation of this software has not been performed. This prototype code is released in open source form with an MIT license (see the LICENSE file in the top level directory of the repository). 

## Getting started

Use the following directions to quickly get started and run ExpressionMatrix2 on a toy dataset distributed as part of the repository.

- Use a Linux machine running Ubuntu 16.04 (the current long term support release of Ubuntu). Other distributions using similar versions of the Linux kernel, such as Linux Mint 18, may also work.

- Make sure the following packages, all available from the standard Ubuntu repositories, are installed:
    * The Boost libraries (package libboost-all-dev).
    * Python3 (package python3-all-dev).
    * Graphviz (package graphviz). 
    * Libraries needed to read HDF5 files (packages libhdf5-10 and libhdf5-cpp-11).

If one or more of these packages are not installed, you can install them using command "apt install packageName", where packageName is as indicated above. This requires root access.

- Select one of the [available releases](https://github.com/chanzuckerberg/ExpressionMatrix2/releases) and download the tar file for the release you selected. 

- Create an empty directory. In the rest of these directions we assume that the directory you created is ~/ExpressionMatrix2. If you use a different name, make sure to use the correct name instead of the one used in the rest of these directions.

- Move the tar file you downloaded from the directory where your browser stores downloads (usually ~/Downloads) to the ~/ExpressionMatrix2 directory your just created, then cd to ~/ExpressionMatrix2.

- Extract the files contained in the tar file you downloaded using command "tar -xvf file.tar", replacing the name of tar file
with the name of the file you downloaded. This will create a doc directory, a test directory, and, most important, a shared library ExpressionMatrix2.so which works as a Python package.

- To make ExpressionMatrix2.so visible as a Python package, set environment variable PYTHONPATH to the name of the directory you are using for your tests, ~/ExpressionMatrix2. 

- Now cd to test directory tests/ToyTest1. It contains csv files with a sample expression matrix and cell meta data for a small toy test with just 3 cells and 3 genes. It also contains two small python scripts. 

- In your copy of the test directory, run command "./run.py". This will read the csv files and create binary data structures stored in mapped files in a new directory named data.

- Run command "./runServer.py". This starts the ExpressionMatrix2 code in a mode where if behaves as an http server that can be used for interactive visualization and analysis. This command continues running until you interrupt it.

- While the server is running, start a web browser on the same machine, and point it to the following URL: http://localhost:17100.

- Interactively explore the functionality offered by the server. You will need a more real test case to perform any interesting visualization or analysis, but if you made it here you are ready to run the ExpressionMatrix2 software on real data. 

You can finds more detailed documentation in the doc directory that was just created by the tar command, ~ExpressionMatrix2/doc. You can point your browser to doc/index.html to get started.

