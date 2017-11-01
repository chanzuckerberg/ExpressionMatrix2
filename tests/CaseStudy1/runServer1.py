#!/usr/bin/python3

from ExpressionMatrix2 import *

# Access our existing expression matrix.
e = ExpressionMatrix(directoryName = 'data')

# Create the graphs we want to see.
# Graphs are not persistent and need to be
# recreated every time the server starts.
e.createCellGraph(
    graphName = 'All', 
    similarPairsName = 'Lsh',
    similarityThreshold = 0.5,
    k=20)
e.createClusterGraph(
    cellGraphName = 'All',
    clusterGraphName = 'All'
    )

# Run the server.
e.explore()



    
    
    
