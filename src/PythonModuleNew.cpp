// This file uses the pybind11 library to export portions of
// the C++ API to Python.
// If something is added or changed here, corresponding documentation
// changes should also be made in doc/PythonApiReference.html.



// Macro that controls exposing to Python of functions declared in filesystem.hpp.
// This is for testing only and should normally be turned off.
#define CZI_EXPRESSION_MATRIX2_TEST_FILESYSTEM 0


// CZI.
#include "ClusterGraph.hpp"
#include "ExpressionMatrix.hpp"
#include "heap.hpp"
#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfLists.hpp"
#include "multipleSetUnion.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#if CZI_EXPRESSION_MATRIX2_TEST_FILESYSTEM
#include "filesystem.hpp"
#endif

// Pybind11.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(ExpressionMatrix2, module)
{
    /*
    module.doc() =
        "Software for analysis, visualization, and clustering "
        "of gene expression data from single-cell RNA sequencing.";
    */

    // Class ExpressionMatrix.
    class_<ExpressionMatrix>(
        module,
        "ExpressionMatrix",
        "This is the top level class in the ExpressionMatrix2 module. "
        "Most high level functionality is provided by this class. "
        "Binary data files for an instance of this class are stored "
        "in a single directory on disk. They are accessed as memory mapped files. ")
       .def(init<string, uint64_t, uint64_t, uint64_t, uint64_t>(),
           "This constructor creates a new (empty) ExpressionMatrix object "
           "in the specified directory. "
           "The directory must not exists. "
           "See `here <../../../Running.html#creationParameters>`__ for more information",
           arg("directoryName"),
           arg("geneCapacity"),
           arg("cellCapacity"),
           arg("cellMetaDataNameCapacity"),
           arg("cellMetaDataValueCapacity")
       )
       .def(init<string, bool>(),
           "This constructor can be used to access an existing ExpressionMatrix object "
           "in the specified directory. The directory must exist. "
           "If write access is not permitted on some of the data, "
           "the operation fails if allowReadOnly is False, "
           "and falls back to read-only access If allowReadOnly is True. "
           "In this case, however, any write operations on data with read-only access "
           "will generate a Segmentation Fault and cause immediate termination "
           "of the calling Python script, without cleanup. "
           "Therefore, this constructor should be called with allowReadOnly = False "
           "except in circumstances where limited functionality "
           "with read-only access to the data is desired. ",
           arg("directoryName"),
           arg("allowReadOnly")=false
       )

       // Get the total number of genes or cells currently in the system.
       .def("geneCount",
           &ExpressionMatrix::geneCount,
           "Returns the total number of genes."
       )
       .def("cellCount",
           &ExpressionMatrix::cellCount,
           "Returns the total number of cells."
       )

       // Genes.
       .def("addGene",
           &ExpressionMatrix::addGene,
           "Adds a gene with the specified name. "
           "Returns True if the gene was added and False if the gene already existed. ",
           arg("geneName")
           )
       .def("geneName",
           &ExpressionMatrix::geneName,
           "Return the name of a gene given its numeric id.",
           arg("geneId")
       )
       .def("geneIdFromName",
           &ExpressionMatrix::geneIdFromName,
           "Return the numeric id of a gene given its name.",
           arg("geneName")
       )


       // Various ways to add cells.
       .def
       (
           "addCell",
           &ExpressionMatrix::addCell,
           "Adds a cell to the system. The cell expression counts "
           "and meta data are given in Python lists. "
           "metadata is a list of tuples with meta data "
           "(name, value) pairs. "
           "expressionCounts is a list of tuples (geneName, count). "
           "See `here <../../../PythonApi.html#addCell>`__ "
           "for an example. "
           "Returns the cell id of the cell that was just added. "
           "Cell ids begin at zero and increment by one each time a cell is added. ",
           arg("metaData"),
           arg("expressionCounts")
       )
       .def
       (
           "addCellFromJson",
           &ExpressionMatrix::addCellFromJson,
           "Adds a cell to the system. The cell expression counts "
           "and meta data are given in a JSON string. "
           "See `here <../../../PythonApi.html#addCellFromJson>`__ "
           "for information on the expected format of the JSON string. "
           "Returns the cell id of the cell that was just added. "
           "Cell ids begin at zero and increment by one each time a cell is added. ",
           arg("jsonString")
       )
       .def("addCells",
           &ExpressionMatrix::addCells,
           "Adds cells to the system, reading expression counts and cell meta data from files. "
           "See `here <../../../PythonApi.html#addCells>`__ for more information. .",
           arg("expressionCountsFileName"),
           arg("expressionCountsFileSeparators") = ",",
           arg("cellMetaDataFileName"),
           arg("cellMetaDataFileSeparators") = ","
       )
#ifndef CZI_EXPRESSION_MATRIX2_SKIP_HDF5
       .def("addCellsFromHdf5",
           &ExpressionMatrix::addCellsFromHdf5,
           "Adds cells to the system, reading expression counts "
           "from a file in HDF5 format as created by the 10X Genomics pipeline. "
           "See `here <../../../PythonApi.html#addCellsFromHdf5>`__ for more information. ",
           arg("fileName"),
           arg("cellNamePrefix"),
           arg("cellMetaData"),
           arg("totalExpressionCountThreshold")
       )
#endif
       .def("addCellsFromBioHub1",
           &ExpressionMatrix::addCellsFromBioHub1,
           "Add to the system cells created by the BioHub pipeline "
           "(July 2017 version, for Illumina data). "
           "See `here <../../../PythonApi.html#addCellsFromBioHub1>`__ for more information. ",
           arg("expressionCountsFileName"),
           arg("initialMetaDataCount"),
           arg("finalMetaDataCount"),
           arg("plateMetaDataFileName")
       )
#ifndef CZI_EXPRESSION_MATRIX2_SKIP_HDF5
       .def("addCellsFromBioHub2",
           &ExpressionMatrix::addCellsFromBioHub2,
           "Add to the system cells created by the BioHub pipeline "
           "(September 2017 version, for 10X data). "
           "See `here <../../../PythonApi.html#addCellsFromBioHub2>`__ for more information. ",
           arg("plateFileName"),
           arg("totalExpressionCountThreshold")
       )
#endif
	   .def("addCellsFromBioHub3",
           &ExpressionMatrix::addCellsFromBioHub3,
           "Add to the system cells created by the BioHub pipeline "
           "(November 2017 version, for Illumina data). "
           "See `here <../../../PythonApi.html#addCellsFromBioHub3>`__ for more information. ",
           arg("expressionCountsFileName"),
           arg("expressionCountsFileSeparators") = ",",
           arg("plateMetaData")
       )
       .def("addCellMetaData",
           &ExpressionMatrix::addCellMetaData,
           "Add cell meta data reading it from a csv file. "
           "Each line of the file contains cell meta data values for one cell. "
           "Cell meta data names are in the first line. "
           "The first column of the file contain cell names. "
           "All cell names must already be present. ",
           arg("cellMetaDataFileName")
       )



       // Accessors for cell meta data.
       .def
       (
           "getCellMetaDataValue",
           (
               string (ExpressionMatrix::*)
               (CellId, const string& metaDataName) const
           )
           &ExpressionMatrix::getCellMetaData,
           "Returns the value of the specified meta data name for the given cell. "
           "Returns an empty string if the cell does not have that meta data name. "
           "The only meta data field that all cells are guaranteed to have is CellName. ",
           arg("cellId"),
           arg("metaDataName")
       )
       .def
       (
           "getCellMetaData",
           (
               vector< pair<string, string> > (ExpressionMatrix::*)
               (CellId) const
           )
           &ExpressionMatrix::getCellMetaData,
           "Returns all meta data pairs (name, value) for a given cell.",
           arg("cellId")
       )
       .def
       (
           "getCellsMetaData",
           (
               vector< vector< pair<string, string> > > (ExpressionMatrix::*)
               (const vector<CellId>&) const
           )
           &ExpressionMatrix::getCellMetaData,
           "Returns all meta data pairs (name, value) "
           "for a set of cells specified by list cellIds. "
           "Each element in the returned list corresponds to the cell id at the same position in list cellIds. ",
           arg("cellIds")
       )
       .def("cellIdFromString",
           &ExpressionMatrix::cellIdFromString,
           "Returns the cell id corresponding to a given name, "
           "or invalidCellId if no cell with the given name exists. ",
           arg("cellString")
       )
       .def("computeMetaDataRandIndex",
           &ExpressionMatrix::computeMetaDataRandIndex,
           "Returns the Rand Index and Adjusted Rand Index "
           "of the contingency table created using two given meta data fields. ",
           arg("cellSetName") = "AllCells",
           arg("metaDataName0"),
           arg("metaDataName1")
       )



       // Accessors for expression counts.
       .def
       (
           "getCellExpressionCount",
           (
               float (ExpressionMatrix::*)
               (CellId, GeneId) const
           )
           &ExpressionMatrix::getCellExpressionCount,
           "Returns the expression count for a given cell and gene. "
           "The returned value can be zero.",
           arg("cellId"),
           arg("geneId")
       )
       .def
       (
           "getCellExpressionCountFromGeneName",
           (
               float (ExpressionMatrix::*)
               (CellId, const string& geneName) const
           )
           &ExpressionMatrix::getCellExpressionCount,
           "Returns the expression count for a given cell and gene. "
           "The returned value can be zero.",
           arg("cellId"),
           arg("geneName")
       )
       .def("getCellExpressionCounts",
           &ExpressionMatrix::getCellExpressionCounts,
           "Returns the non-zero expression counts for a given cell. "
           "The returned list contains pairs of gene ids "
           "and the corresponding expression counts for the given cell.",
           arg("cellId")
       )
       .def
       (
           "getCellsExpressionCount",
           (
               vector<float> (ExpressionMatrix::*)
               (const vector<CellId>&, GeneId) const
           )
           &ExpressionMatrix::getCellsExpressionCount,
           "Returns the expression counts for a set of cells and for a given gene. "
           "Some or all of the returned values can be zero. ",
           arg("cellIds"),
           arg("geneId")
       )
       .def
       (
           "getCellsExpressionCountFromGeneName",
           (
               vector<float> (ExpressionMatrix::*)
               (const vector<CellId>&, const string& geneName) const
           )
           &ExpressionMatrix::getCellsExpressionCount,
           "Returns the expression counts for a set of cells and for a given gene. "
           "Some or all of the returned values can be zero. ",
           arg("cellIds"),
           arg("geneName")
       )
       .def("getCellsExpressionCounts",
           &ExpressionMatrix::getCellsExpressionCounts,
           "Returns the non-zero expression counts for a given set of cells. "
           "The returned list contains, for each cell, "
           "a list of pairs of gene ids and the corresponding expression counts.",
           arg("cellIds")
           )
       .def("getCellsExpressionCountsForGenes",
           &ExpressionMatrix::getCellsExpressionCountsForGenes,
           "Returns the non-zero expression counts for given sets of cells and genes. "
           "The returned list contains, for each cell, pairs of gene ids "
           "and the corresponding expression counts.",
           arg("cellIds"),
           arg("geneIds")
       )



       // Gene sets.
       .def(
           "createGeneSetUsingInformationContent",
           (
               void (ExpressionMatrix::*)
               (
                   const string& existingGeneSetName,
                   const string& cellSetName,
                   NormalizationMethod normalizationMethod,
                   double geneInformationContentThreshold,
                   const string& newGeneSetName
               )
           )
           &ExpressionMatrix::createGeneSetUsingInformationContent,
           "Creates a new gene set consisting of genes that pass an information content test. "
           "The information content, in bits, of each gene in existingGeneSetName "
           "is computed using cellSetName and the specified NormalizationMethod. "
           "The new gene set newGeneSetName is created, "
           "consisting of genes whose computed information content (in bits) "
           "exceeds the specified geneInformationContentThreshold. "
           "The information content is roughly related to the variability "
           "of the expression of that gene in the specified cell set. "
           "Genes that are widely expressed have a low information content. "
           "Genes that are very selectively expressed have a high information content. "
           "A gene that is equally expressed in all cells has zero information content. "
           "A gene that is expressed only in one cell has log2N bits of information content, "
           "where N is the number of cells in the cell set used to compute information content. "
           "A geneInformationContentThreshold of 2 bits is often effective "
           "in filtering out widely expressed genes. ",
           arg("existingGeneSetName") = "AllGenes",
           arg("cellSetName") = "AllCells",
           arg("normalizationMethod"),
           arg("geneInformationContentThreshold"),
           arg("newGeneSetName")
       )
       .def("createGeneSetIntersection", &ExpressionMatrix::createGeneSetIntersection,
           "Creates a new gene set as the intersection of two or more gene sets. "
           "The names of the input gene sets specified as the first argument "
           "must be separated by commas. "
           "Returns False if the operation failed because one of the input gene sets "
           "does not exist or the output gene set already exists, True otherwise. ",
           arg("inputSetsNames"),
           arg("outputSetName")
       )
       .def("createGeneSetUnion", &ExpressionMatrix::createGeneSetUnion,
           "Creates a new gene set as the union of two or more gene sets. "
           "The names of the input gene sets specified as the first argument "
           "must be separated by commas. "
           "Returns False if the operation failed because one of the input gene sets "
           "does not exist or the output gene set already exists, True otherwise.",
           arg("inputSetsNames"),
           arg("outputSetName")
       )
       .def
       (
           "createGeneSetDifference",
           (
               bool (ExpressionMatrix::*)
               (const string&, const string&, const string&)
           )
           &ExpressionMatrix::createGeneSetDifference,
           "Creates a new gene set as the difference between two gene sets. "
           "The newly created gene set consists of all genes that are in inputGeneSetName0 "
           "but not in inputGeneSetName1. "
           "Returns False if the operation failed because one of the input gene sets "
           "does not exist or the output gene set already exists, True otherwise. ",
           arg("inputSetName0"),
           arg("inputSetName1"),
           arg("oututSetName")
       )
       .def
       (
           "createGeneSetFromGeneNames",
           (
               void (ExpressionMatrix::*)
               (const string&, const vector<string>&)
           )
           &ExpressionMatrix::createGeneSetFromGeneNames,
           "Creates a new gene set given a list of gene names.",
           arg("geneSetName"),
           arg("geneNames")
       )



        // Cell sets.
       .def
       (
           "createCellSetUsingMetaData",
           (
               bool (ExpressionMatrix::*)
               (const string&, const string&, const string&, bool)
           )
           &ExpressionMatrix::createCellSetUsingMetaData,
           "Creates a new cell set cellSetName consisting of all cells "
           "for which the specified metaDataName matches the given matchTarget. "
           "If useRegularExpression is False, this uses a simple string matching. "
           "Otherwise, matchTarget is treated as a regular expression to be matched. "
           "As a simple application example, suppose all cells contain a meta data name Tissue "
           "that gives the type of tissue from which the cell was obtained. "
           "Calling createCellSetUsingMetaData with metaDataName='Tissue', matchTarget='Brain', "
           "and useRegularExpression='False' will create a new cell set "
           "consisting of the cells for which the value of the Tissue meta data field is Brain. "
           "Regular expressions can be used to select values that satisfy a variety of criteria.",
           arg("cellSetName"),
           arg("metaDataFieldName"),
           arg("matchString"),
           arg("useRegex")
       )
       .def("createCellSetUsingNumericMetaDataGreaterThan",
           &ExpressionMatrix::createCellSetUsingNumericMetaDataGreaterThan,
           "Creates a new cell set cellSetName consisting of all cells "
           "for which the specified metaDataName is numeric and greater than lowerBound.",
           arg("cellSetName"),
           arg("metaDataFieldName"),
           arg("lowerBound")
       )
       .def("createCellSetUsingNumericMetaDataLessThan",
           &ExpressionMatrix::createCellSetUsingNumericMetaDataLessThan,
           "Creates a new cell set cellSetName consisting of all cells "
           "for which the specified metaDataName is numeric and less than upperBound. "
           "",
           arg("cellSetName"),
           arg("metaDataFieldName"),
           arg("upperBound")
       )
       .def("createCellSetUsingNumericMetaDataBetween",
           &ExpressionMatrix::createCellSetUsingNumericMetaDataBetween,
           "Creates a new cell set cellSetName consisting of all cells "
           "for which the specified metaDataName is numeric, "
           "greater than lowerBound, and less than upperBound. ",
           arg("cellSetName"),
           arg("metaDataFieldName"),
           arg("lowerBound"),
           arg("upperBound")
       )
       .def("createCellSetIntersection",
           &ExpressionMatrix::createCellSetIntersection,
           "Creates a new cell set as the intersection of two or more cell sets. "
           "The names of the input cell sets specified as the first argument "
           "must be separated by commas. "
           "Returns False if the operation failed because one of the input cell sets does not exist "
           "or the output cell set already exists, True otherwise. ",
           arg("inputSetsNames"),
           arg("outputSetName")
       )
       .def("createCellSetUnion",
           &ExpressionMatrix::createCellSetUnion,
           "Creates a new cell set as the union of two or more cell sets. "
           "The names of the input cell sets specified as the first argument "
           "must be separated by commas. "
           "Returns False if the operation failed because one of the input cell sets does not exist "
           "or the output cell set already exists, True otherwise. ",
           arg("inputSetsNames"),
           arg("outputSetName")
       )
       .def
       (
           "createCellSetDifference",
           (
               bool (ExpressionMatrix::*)
               (const string&, const string&, const string&)
           )
           &ExpressionMatrix::createCellSetDifference,
           "Creates a new cell set as the difference between two cell sets. "
           "The newly created cell set consists of all cells "
           "that are in inputCellSetName0 but not in inputCellSetName1. "
           "Returns False if the operation failed because one of the input cell sets does not exist "
           "or the output cell set already exists, True otherwise. ",
           arg("inputSetName0"),
           arg("inputSetName1"),
           arg("outputSetName")
       )
       .def("getCellSetNames",
           &ExpressionMatrix::getCellSetNames,
           "Returns a list containing the names of all currently defined cell sets. "
       )
       .def("getCellSet",
           &ExpressionMatrix::getCellSet,
           "Returns the cell ids of the cells in the cell set with the specified name. "
           "If the cell set does not exists, returns an empty container.",
           arg("cellSetName")
       )
       .def
       (
           "removeCellSet",
           (
               void (ExpressionMatrix::*)
               (const string&)
           )
           &ExpressionMatrix::removeCellSet,
           "Removes the cell set with the specified name. "
           "If the cell set does not exists, raises an exception.",
           arg("cellSetName")
       )



       // Compute cell similarity.
       .def("computeCellSimilarity",
           (
               double (ExpressionMatrix::*)
               (const string&, CellId, CellId) const
           )
           &ExpressionMatrix::computeCellSimilarity,
           "Returns the similarity between the expression vectors of two cells, "
           "specified by their cell ids. "
           "The similarity is computed taking into account only genes in the specifified gene set. "
           "It is computed as the `Pearson correlation coefficient "
           "<https://en.wikipedia.org/wiki/Pearson_correlation_coefficient>`__ "
           "of the expression vectors of the two cells. "
           "The computed similarity is always between -1 and 1, and it is rarely negative. ",
           arg("geneSetName") = "AllGenes",
           arg("cellId0"),
           arg("cellId1")
       )
       .def("computeApproximateLshCellSimilarity",
           &ExpressionMatrix::computeApproximateLshCellSimilarity,
           "Only intended to be used for testing. "
           "See the source code in the ExpressionMatrix2/src directory for more information. "
       )
       .def("writeLshSimilarityComparison",
           &ExpressionMatrix::writeLshSimilarityComparison,
           "Only intended to be used for testing. "
           "See the source code in the ExpressionMatrix2/src directory for more information. "
       )
       .def("writeLshSimilarityComparisonSlow",
           &ExpressionMatrix::writeLshSimilarityComparisonSlow,
           "Only intended to be used for testing. "
           "See the source code in the ExpressionMatrix2/src directory for more information. "
       )



       // Find similar pairs of cells.
       .def("findSimilarPairs0",
           &ExpressionMatrix::findSimilarPairs0,
           "Creates and stores a new object similarPairsName "
           "containing pairs of similar cells from cellSetName. "
           "The computation is done by directly looping "
           "over all possible cell pairs with both cells in cellSetName "
           "and performing for each pair an exact computation of the cell similarity, "
           "taking into account only genes in geneSetName. "
           "For each cell, only the best (most similar) k or fewer similar cells are stored, "
           "among those that exceed the specified similarityThreshold. "
           "This computational cost of this function grows "
           "with the square of the number of cells in cellSetName. "
           "The required computing time will typically be a few minutes "
           "for a few thousand cells or several hours for a few tens of thousands cells. "
           "When the number of cells exceeds a few thousands, "
           "it is more practical to perform an approximate computation using findSimilarPairs4. ",
           arg("geneSetName") = "AllGenes",
           arg("cellSetName") = "AllCells",
           arg("similarPairsName"),
           arg("k") = 100,
           arg("similarityThreshold") = 0.2
       )
       .def("findSimilarPairs1",
           &ExpressionMatrix::findSimilarPairs1,
           "Obsolete. Use findSimilarPairs4 instead."
       )
       .def("findSimilarPairs2",
           &ExpressionMatrix::findSimilarPairs2,
           "Obsolete. Replaced by findSimilarPairs5."
           "It will provide LSH-based computation of similar cell pairs "
           "without looping over all possible pairs of cells."
       )
       .def("findSimilarPairs3",
           &ExpressionMatrix::findSimilarPairs3,
           "Obsolete. Use findSimilarPairs4 instead.",
           arg("geneSetName") = "AllGenes",
           arg("cellSetName") = "AllCells",
           arg("similarPairsName"),
           arg("k") = 100,
           arg("similarityThreshold") = 0.2,
           arg("lshCount") = 1024,
           arg("seed") = 231
       )
       .def("findSimilarPairs3Benchmark",
           &ExpressionMatrix::findSimilarPairs3Benchmark,
           "Like findSimilarPairs3, but does not store anything. Used for benchmarks.",
           arg("geneSetName") = "AllGenes",
           arg("cellSetName") = "AllCells",
           arg("lshCount") = 1024,
           arg("seed") = 231
       )
       .def("findSimilarPairs4",
           &ExpressionMatrix::findSimilarPairs4,
           "Like findSimilarPairs0, but uses Locality-Sensitive Hashing (LSH) "
           "with lshCount LSH hashes to speed up the computation of the similar pairs. "
           "Like findSimilarPairs0, this still loops over all possible cell pairs "
           "with both cells in cellSetName and therefore the computational cost "
           "still grows asymptotically with the square of the number of cells. "
           "However the computation is orders of magnitudes faster. "
           "The computation is approximate, and the error decreases as lshCount increases. "
           "For the suggested value lshCount=1024, "
           "the standard deviation of the pair similarity computed in this way is 0.05 or less.",
           arg("geneSetName") = "AllGenes",
           arg("cellSetName") = "AllCells",
           arg("similarPairsName"),
           arg("k") = 100,
           arg("similarityThreshold") = 0.2,
           arg("lshCount") = 1024,
           arg("seed") = 231
       )
       .def("findSimilarPairs5",
           &ExpressionMatrix::findSimilarPairs5,
           "LSH-based computation of similar cell pairs "
           "without looping over all possible pairs of cells.",
           arg("geneSetName") = "AllGenes",
           arg("cellSetName") = "AllCells",
           arg("similarPairsName"),
           arg("k") = 100,
           arg("similarityThreshold") = 0.2,
           arg("lshCount") = 1024,
           arg("lshSliceLength"),
           arg("bucketOverflow") = 1000,
           arg("seed") = 231
       )
       .def("writeSimilarPairs",
           &ExpressionMatrix::writeSimilarPairs,
           "Only intended to be used for testing. "
           "See the source code in the ExpressionMatrix2/src directory for more information. "
           "For debugging/testing use only."
       )
       .def("analyzeSimilarPairs",
           &ExpressionMatrix::analyzeSimilarPairs,
           "Only intended to be used for testing. "
           "See the source code in the ExpressionMatrix2/src directory for more information. "
       )
       .def("compareSimilarPairs",
           &ExpressionMatrix::compareSimilarPairs,
           "Only intended to be used for testing. "
           "See the source code in the ExpressionMatrix2/src directory for more information. "
       )
       .def("analyzeLsh",
           &ExpressionMatrix::analyzeLsh,
           "Only intended to be used for testing. "
           "See the source code in the ExpressionMatrix2/src directory for more information. "
       )
       .def("analyzeLshSignatures",
           &ExpressionMatrix::analyzeLshSignatures,
           "Only intended to be used for testing. "
           "See the source code in the ExpressionMatrix2/src directory for more information. ",
           arg("geneSetName") = "AllGenes",
           arg("cellSetName") = "AllCells",
           arg("lshCount") = 1024,
           arg("seed") = 231
       )



       // Cell graphs.
       .def("getCellGraphNames",
           &ExpressionMatrix::getCellGraphNames,
           "Return a list containing the names of all currently defined cell similarity graphs."
       )
       .def
       (
           "createCellGraph",
           (
               void (ExpressionMatrix::*)
               (const string&, const string&, const string&, double, size_t)
           )
           &ExpressionMatrix::createCellGraph,
           "Create a new cell similarity graph. "
           "similarityThreshold is the minimum cell similarity required for an edge to be created. "
           "k controls the maximum connectivity of the graph. "
           "A cell graph or cell similarity graph is an undirected graph "
           "in which each vertex represents a cell. "
           "An edge between two vertices is created if the corresponding cells "
           "have similarity that exceeds a chosen threshold. "
           "In areas of high connectivity, only the edges "
           "with the highest similarity for each cell are kept. ",
           arg("graphName"),
           arg("cellSetName") = "AllCells",
           arg("similarPairsName"),
           arg("similarityThreshold") = 0.5,
           arg("k") = 20
       )
       .def("computeCellGraphLayout",
           &ExpressionMatrix::computeCellGraphLayout,
           "Computes the two-dimensional layout for the graph with the given name. "
           "The graph must have been previously created with a call to "
           ":py:func:`ExpressionMatrix2.ExpressionMatrix.createCellGraph`. "
           "This function must be called once before the first call to "
           ":py:func:`ExpressionMatrix2.ExpressionMatrix.getCellGraphVertices` for the graph. ",
           arg("graphName")
       )
       .def("getCellGraphVertices",
           &ExpressionMatrix::getCellGraphVertices,
           "Returns information about the vertices of the cell graph with the given name. "
           "The graph must have been previously created with a call to "
           ":py:func:`ExpressionMatrix2.ExpressionMatrix.createCellGraph` "
           "and its two-dimensional layout must have been computed with a call to "
           ":py:func:`ExpressionMatrix2.ExpressionMatrix.computeCellGraphLayout`. ",
           arg("graphName")
       )
       .def("getCellGraphEdges",
           &ExpressionMatrix::getCellGraphEdges,
           "Returns the cell ids corresponding to the two vertices "
           "of each edge of the cell graph with the given name. "
           "The graph must have been previously created with a call to "
           ":py:func:`ExpressionMatrix2.ExpressionMatrix.createCellGraph`.",
           arg("graphName")
       )



       // Cluster graphs.
       .def
       (
           "createClusterGraph",
           (
               void (ExpressionMatrix::*)
               (const string&, const string&, size_t, size_t, size_t, size_t, size_t, double)
           )
           &ExpressionMatrix::createClusterGraph,
           "Creates a new cluster graph by running label propagation clustering "
           "on an existing cell graph. "
           "A cluster graph is an undirected graph "
           "in which each vertex represents a cluster found in a cell graph. "
           "An edge between two vertices is created if the corresponding clusters "
           "have average expression vectors that are sufficiently similar. "
           "In areas of high connectivity, only the edges "
           "with the highest similarity for each cluster are kept. ",
           arg("cellGraphName"),
           arg("clusterGraphName"),
           arg("stableIterationCount") = 3,
           arg("maxIterationCount") = 100,
           arg("seed") = 231,
           arg("minClusterSize") = 100,
           arg("k") = 3,
           arg("similarityThreshold") = 0.5
       )
       .def("getClusterGraphVertices",
           &ExpressionMatrix::getClusterGraphVertices,
           "Returns a list of the cluster ids for the vertices of an existing cluster graph. "
           "Each vertex of a cluster graph corresponds to a cluster. ",
           arg("clusterGraphName")
       )
       .def("getClusterGraphGenes",
           &ExpressionMatrix::getClusterGraphGenes,
           "Returns a list of the ids of the genes used to create a cluster graph. ",
           arg("clusterGraphName")
       )
       .def("getClusterCells",
           &ExpressionMatrix::getClusterCells,
           "Returns a list of the cell ids for the cells a cluster. "
           "The cluster is specified using the name of the cluster graph it belongs to, "
           "and its cluster id within that cluster graph. ",
           arg("clusterGraphName"),
           arg("clusterId")
       )
       .def("getClusterAverageExpression",
           &ExpressionMatrix::getClusterAverageExpression,
           "Return the average, L2-normalized, gene expression of a cluster. "
           "The cluster is specified using the name of the cluster graph it belongs to, "
           "and its cluster id within that cluster graph. "
           "The values returned correspond one on one to the gene ids "
           "returned by :py:func:`ExpressionMatrix2.ExpressionMatrix.getClusterGraphGenes` for the same cluster graph. .",
           arg("clusterGraphName"),
           arg("clusterId")
       )
       .def
       (
           "createMetaDataFromClusterGraph",
           (
               void (ExpressionMatrix::*)
               (const string&, const string&)
           )
           &ExpressionMatrix::createMetaDataFromClusterGraph,
           "Creates meta data from the cluster ids stored in a cluster graph by "
           "storing the cluster ids in the specified cluster graph into a cell meta data field.",
           arg("clusterGraphName"),
           arg("metaDataName")
       )



       // Run the http server.
       .def("explore",
           (
               void (ExpressionMatrix::*)
               (uint16_t, const string&)
           )
           &ExpressionMatrix::explore,
           "Starts an http server that can be used, in conjunction with a Web browser, "
           "to interact with the ExpressionMatrix object. ",
           arg("port") = 17100,
           arg("docDirectory") = ""
       )



       // Functions used only for testing or debugging.
       .def("testExpressionMatrixSubset",
           &ExpressionMatrix::testExpressionMatrixSubset,
           "For debugging/testing use only."
       )
       ;



    // Class CellGraphVertexInfo.
    class_<CellGraphVertexInfo>(
        module,
        "CellGraphVertexInfo",
        "This class is used to give the Python code access to information "
        "about the vertices of a cell similarity graph. "
        "An object of this type contains information for a single vertex.")
        .def(init<>(), "This constructor creates an uninitialized CellGraphVertexInfo object.")
        .def_readonly("cellId", &CellGraphVertexInfo::cellId,
            "This data member contains the cell id of the cell corresponding to the vertex. ")
        .def("x", &CellGraphVertexInfo::x,
            "Returns the x coordinate of the vertex in the two-dimensional layout "
            "of the graph the vertex belongs to. ")
        .def("y", &CellGraphVertexInfo::y,
            "Returns the y coordinate of the vertex in the two-dimensional layout "
            "of the graph the vertex belongs to. ")
        ;



    // Enum class NormalizationMethod.
    enum_<NormalizationMethod>(
        module,
        "NormalizationMethod",
        "Various ways to normalize gene expressions.\n\n"
        "- none: no normalization.\n"
        "- L1: L1 normalization (sum of values is 1).\n"
        "- L2: L2 normalization (sum of squares of values is 1).\n"
        "- Invalid: invalid normalization.\n"
        )
        .value(normalizationMethodToShortString(NormalizationMethod::None).c_str(),       NormalizationMethod::None)
        .value(normalizationMethodToShortString(NormalizationMethod::L1).c_str(),         NormalizationMethod::L1)
        .value(normalizationMethodToShortString(NormalizationMethod::L2).c_str(),         NormalizationMethod::L2)
        .value(normalizationMethodToShortString(NormalizationMethod::Invalid).c_str(),    NormalizationMethod::Invalid)
        .export_values()
        ;


    // Some non-member functions used only for testing or debugging.
    module.def("testMemoryMappedVector",
        testMemoryMappedVector,
        "Only intended to be used for testing. "
        "See the source code in the ExpressionMatrix2/src directory for more information. "
        );
    module.def("testMemoryMappedVectorOfLists",
        testMemoryMappedVectorOfLists,
        "Only intended to be used for testing. "
        "See the source code in the ExpressionMatrix2/src directory for more information. "
        );
    module.def("testMemoryMappedStringTable",
        testMemoryMappedStringTable,
        "Only intended to be used for testing. "
        "See the source code in the ExpressionMatrix2/src directory for more information. "
        );
    module.def("testHeap",
        testHeap,
        "Only intended to be used for testing. "
        "See the source code in the ExpressionMatrix2/src directory for more information. "
        );
    module.def("testKeepBest",
        testKeepBest,
        "Only intended to be used for testing. "
        "See the source code in the ExpressionMatrix2/src directory for more information. "
        );
    module.def("multipleSetUnionTest",
        multipleSetUnionTest,
        "Only intended to be used for testing. "
        "See the source code in the ExpressionMatrix2/src directory for more information. "
        );



#if CZI_EXPRESSION_MATRIX2_TEST_FILESYSTEM
    // Functions declared in filesystem.hpp.
    // This is for testing only and should normally be turned off.
    module.def("exists", ExpressionMatrix2::filesystem::exists);
    module.def("isRegularFile", ExpressionMatrix2::filesystem::isRegularFile);
    module.def("isDirectory", ExpressionMatrix2::filesystem::isDirectory);
    module.def("createDirectory", ExpressionMatrix2::filesystem::createDirectory);
    module.def("remove", ExpressionMatrix2::filesystem::remove);
    module.def("directoryContents", ExpressionMatrix2::filesystem::directoryContents);
    module.def("extension", ExpressionMatrix2::filesystem::extension);
    module.def("fileName", ExpressionMatrix2::filesystem::fileName);
#endif
}
