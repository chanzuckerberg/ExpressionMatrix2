#!/usr/bin/python3



import ExpressionMatrix2 


# Construct the ExpressionMatrix using the binary data stored in directory "data".
e = ExpressionMatrix2.ExpressionMatrix('data')

# Run the server.
serverParameters = ExpressionMatrix2.ServerParameters()
serverParameters.port = 17100
e.explore(serverParameters)



    
    
    
