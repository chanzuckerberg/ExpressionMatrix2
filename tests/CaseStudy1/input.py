#!/usr/bin/python3


from ExpressionMatrix2 import *

# Create a new, empty expression matrix.
# The data directory must not exist.
e = ExpressionMatrix(
    directoryName = 'data', 
    geneCapacity = 100000,
    cellCapacity = 10000,
    cellMetaDataNameCapacity = 10000,
    cellMetaDataValueCapacity = 1000000
    )

# Add the cells.
e.addCells(
    expressionCountsFileName = 'GBM_raw_gene_counts.csv',
    expressionCountsFileSeparators = ' ', 
    cellMetaDataFileName = 'GBM_metadata.csv',
    cellMetaDataFileSeparators = ' ' 
    )

print('Input completed.')

 
