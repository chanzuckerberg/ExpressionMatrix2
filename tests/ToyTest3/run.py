#!/usr/bin/python3


# Import the shared library, which behaves as a Python module.
import ExpressionMatrix2 
import json



# Create the expression matrix and add the cells.
# This creates directory "data" to contain the binary data for this expression matrix.
# Later, we can access the binary data using a different ExpressinMatrix constructor (see runServer.py)
e = ExpressionMatrix2.ExpressionMatrix(directoryName = 'data')

cell0 = {'metaData': {'CellName': 'cell0', 'Type': 'type0'}, 'expressionCounts': {'gene1': 10,'gene2': 20}}
cell1 = {'metaData': {'CellName': 'cell1', 'Type': 'type1'}, 'expressionCounts': {'gene1': 30}}
cell2 = {'metaData': {'CellName': 'cell2', 'Type': 'type2'}, 'expressionCounts': {'gene2': 40}}

e.addCellFromJson(jsonString = json.dumps(cell0)) 
e.addCellFromJson(jsonString = json.dumps(cell1)) 
cell2Metadata=[('CellName', 'cell2'), ('Type', 'type2')]
cell2ExpressionCounts = [('gene2', 40.)]
e.addCell(metaData=cell2Metadata, expressionCounts=cell2ExpressionCounts) 
print('There are %i genes and %i cells.' % (e.geneCount(), e.cellCount()))



# Find pairs of similar cells.
e.findSimilarPairs0(similarPairsName = 'Exact',)




    
    
    
