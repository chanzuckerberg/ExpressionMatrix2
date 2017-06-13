# ExpressionMatrix2
Software for analysis of gene expression data from single-cell RNA sequencing.

This repository contains software for analysis of gene expression data from single-cell RNA sequencing developed at Chan-Zuckerberg Initiative.

This code is at a prototype, pre-alpha stage. It is hoped that it can be useful in its current form, but it is likely that it contains bugs or other errors that impact its ability to give meaningful results. In addition, a validation of this code has not been performed. This prototype code is released in open source form with an MIT license (see LICENSE.md). 

The next section contains some information to get you started quickly. For more complete documentation, make a local clone of the repository and point your browser to doc/index.html. An easy way to make a local copy of the repositor is to use one of the click on the green "Clone or download" button for this GitHub repository and choose one of the options available there.

## Getting started

Use the following direction to get quickly started and run ExpressionMatrix2 on a toy dataset distributed as part of the repository.

    Use a Linux machine running Ubuntu 16.04 (the current long term support release of Ubuntu). Other distributions using similar versions of the Linux kernel, such as Linux Mint 18, will also work.
    Make sure the following packages, all available from the standard Ubuntu repositories, are installed:
        - The Boost libraries (package libboost-all-dev).
        - Python3 (package python3-all-dev).
        - Graphviz (package graphviz). 
    If one or more of these packages are not installed, you can install them using command "apt install packageName", where packageName is as indicated above. This requires root access.
    Create an empty directory for the test you are going to run.
    In GitHub, navigate to the bin directory of the chanzuckerberg/ExpressionMatrix2 repository and down to the library you want to use. There may be more than one version, each named ExpressionMatrix2.so and located in a directory named from the date the library was generated. Click on the selected library, then "Download". Copy it to the empty directory you created for the test. (If you prefer, you can instead clone the repository and get the library from you cloned copy. This would require running command "git clone https://github.com/chanzuckerberg/ExpressionMatrix2").
    In the same way, copy to your test directory files ExpressionMatrix.csv, MetaData.csv, run.py, and runServer.py from the GitHub repository under tests/ToyTest1. These files contain the expression matrix and cell meta data for the test.
    Run command "./run.py". This will read the csv files and create binary data structures stored in mapped files in a new directory named data.
    Run command "./runServer.py". This starts the ExpressionMatrix2 code in a mode where if behaves as an http server that can be used for interactive visualization and analysis. This command continues running until you interrupt it.
    While the server is running, start a web browser on the same machine, and point it to the following URL: http://localhost:17100.
    Interactively explore the functionality offered by the server. You will need a more real test case to perform any interesting visualization or analysis, but if you made it here you are ready to run the ExpressionMatrix2 on real data. 

Note that to simplify the "Getting started" process we copied ExpressionMatrix2.so to the same directory containing the Python scripts, so it is found automatically. In real use, this will not be practical, and instead you will want to set environmnet variable PYTHONPATH to the directory containing ExpressionMatrix2.so. This will allow the Python interpreter to locate it. 

