#!/usr/bin/python3


# Import the shared library, which behaves as a Python module.
import ExpressionMatrix2 



# Create the expression matrix and add the cells.
# This creates directory "data" to contain the binary data for this expression matrix.
# Later, we can access the binary data using a different ExpressinMatrix constructor (see runServer.py)
e = ExpressionMatrix2.ExpressionMatrix(
    directoryName = 'data',
    geneCapacity = 1<<18,                # Maximum number of genes.
    cellCapacity = 1<<16,                # Maximum number of cells.           
    cellMetaDataNameCapacity = 1<<12,    # Maximum number of distinct cell meta data name strings.
    cellMetaDataValueCapacity = 1<<20    # Maximum number of distinct cell meta data value strings.
    )
e.addCells(
    expressionCountsFileName = 'ExpressionMatrix1.csv', 
    cellMetaDataFileName = 'MetaData1.csv'
    )
e.addCells(
    expressionCountsFileName = 'ExpressionMatrix2.csv', 
    cellMetaDataFileName = 'MetaData2.csv'
    )



# Find pairs of similar cells.
e.findSimilarPairs0(similarPairsName = 'Exact')




    
    
    
