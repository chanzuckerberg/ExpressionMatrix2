#!/usr/bin/python3


# Import the shared library, which behaves as a Python module.
import ExpressionMatrix2 



# Create the expression matrix and add the cells.
# This creates directory "data" to contain the binary data for this expression matrix.
# Later, we can access the binary data using a different ExpressinMatrix constructor (see runServer.py)
e = ExpressionMatrix2.ExpressionMatrix(directoryName = 'data')
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




    
    
    
