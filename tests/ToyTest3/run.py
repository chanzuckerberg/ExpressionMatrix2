#!/usr/bin/python3


# Import the shared library, which behaves as a Python module.
import ExpressionMatrix2 
import json



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

cell0 = {'metaData': {'CellName': 'cell0', 'Type': 'type0'}, 'expressionCounts': {'gene1': 10,'gene2': 20}}
cell1 = {'metaData': {'CellName': 'cell1', 'Type': 'type1'}, 'expressionCounts': {'gene1': 30}}
cell2 = {'metaData': {'CellName': 'cell2', 'Type': 'type2'}, 'expressionCounts': {'gene2': 40}}

e.addCell(jsonString = json.dumps(cell0)) 
e.addCell(jsonString = json.dumps(cell1)) 
e.addCell(jsonString = json.dumps(cell2)) 
print('There are %i genes and %i cells.' % (e.geneCount(), e.cellCount()))



# Find pairs of similar cells.
e.findSimilarPairs0(similarPairsName = 'Exact',)




    
    
    
