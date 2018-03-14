#!/usr/bin/python3


from ExpressionMatrix2 import *

# Create a new, empty expression matrix.
# The data directory must not exist.
e = ExpressionMatrix(directoryName = 'data')

# Add the cells.
e.addCells(
    expressionCountsFileName = 'GBM_raw_gene_counts.csv',
    expressionCountsFileSeparators = ' ', 
    cellMetaDataFileName = 'GBM_metadata.csv',
    cellMetaDataFileSeparators = ' ' 
    )

print('Input completed.')

 
