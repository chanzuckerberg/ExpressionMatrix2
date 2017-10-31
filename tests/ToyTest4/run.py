#!/usr/bin/python3


# Import the shared library, which behaves as a Python module.
import ExpressionMatrix2 



# Some expression matrix creation parameters.
# The capacity values defined below are hard limits.
# For good hash table performance, the actual numbers of items
# used should stay below half the values specified here.
parameters = ExpressionMatrix2.ExpressionMatrixCreationParameters()
parameters.geneCapacity = 1<<18                 # Maximum number of genes.
parameters.cellCapacity = 1<<16                 # Maximum number of cells.           
parameters.cellMetaDataNameCapacity = 1<<12     # Maximum number of distinct cell meta data name strings.
parameters.cellMetaDataValueCapacity = 1<<20    # Maximum number of distinct cell meta data value strings.



# Create the expression matrix and add the cells.
# This creates directory "data" to contain the binary data for this expression matrix.
# Later, we can access the binary data using a different ExpressinMatrix constructor (see runServer.py)
e = ExpressionMatrix2.ExpressionMatrix(
    directoryName = 'data', 
    parameters = parameters)
e.addCells(
    expressionCountsFileName = 'ExpressionMatrix.csv', 
    cellMetaDataFileName = 'MetaData.csv'
    )



# Find pairs of similar cells.
e.findSimilarPairs0(similarPairsName = 'Exact')




    
    
    
