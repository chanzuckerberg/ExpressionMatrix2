.. ExpressionMatrix2 documentation master file, created by
   sphinx-quickstart on Fri Oct 27 12:48:09 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ExpressionMatrix2 Python API
============================

Contents:

.. toctree::
   :maxdepth: 2



.. automodule:: ExpressionMatrix2



CellGraphVertexInfo
-------------------------

.. autoclass:: CellGraphVertexInfo
   :members: 
   
   .. automethod:: __init__
   .. autoattribute:: cellId



ClusterGraphCreationParameters
------------------------------------

.. autoclass:: ClusterGraphCreationParameters
   :members: 
   
   .. automethod:: __init__
   .. autoattribute:: stableIterationCount 
   .. autoattribute:: maxIterationCount  
   .. autoattribute:: seed   
   .. autoattribute:: minClusterSize    
   .. autoattribute:: maxConnectivity    
   .. autoattribute:: similarityThreshold    



ExpressionMatrix
----------------------

.. autoclass:: ExpressionMatrix
   :members: 
   
   .. automethod:: __init__

.. autoclass:: ExpressionMatrixCreationParameters
   :members: 
   
   .. automethod:: __init__
   .. autoattribute:: geneCapacity  
   .. autoattribute:: cellCapacity   
   .. autoattribute:: cellMetaDataNameCapacity    
   .. autoattribute:: cellMetaDataValueCapacity     



NormalizationMethod
-------------------------

.. autoclass:: NormalizationMethod
   :members: 



ServerParameters
----------------------

.. autoclass:: ServerParameters
   :members: 
   
   .. automethod:: __init__
   .. autoattribute:: port  
   .. autoattribute:: docDirectory    



Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

