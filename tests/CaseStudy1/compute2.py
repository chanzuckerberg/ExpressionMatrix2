#!/usr/bin/python3


from ExpressionMatrix2 import *

# Access our existing expression matrix.
e = ExpressionMatrix(directoryName = 'data')

# Find pairs of similar cell (approximate computation using LSH, 
# while still looping over all pairs). 
print("Finding pairs of similar cells.")
e.findSimilarPairs3(geneSetName='HighInformationGenes', similarPairsName = 'LshHighInfo')





    
    
