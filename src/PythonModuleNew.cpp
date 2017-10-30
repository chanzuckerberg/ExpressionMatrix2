// This file uses the pybind11 library to export portions of
// the C++ API to Python.
// If something is added or changed here, corresponding documentation
// changes should also be made in doc/PythonApiReference.html.


// CZI.
#include "ClusterGraph.hpp"
#include "ExpressionMatrix.hpp"
#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfLists.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

// Pybind11.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(ExpressionMatrix2, module)
{
    module.doc() =
        "Software for analysis, visualization, and clustering "
        "of gene expression data from single-cell RNA sequencing.";

    // Class ExpressionMatrix.
    class_<ExpressionMatrix>(
        module,
        "ExpressionMatrix",
        "Top level object used to store and manipulate an expression matrix "
        "containing information on its genes and cells "
        "and many related data structures.")
       .def(init<string, ExpressionMatrixCreationParameters>(),
           "Construct a new, empty ExpressionMatrix.",
           arg("directoryName"),
           arg("parameters")
       )
       .def(init<string, bool>(),
           "Access an existing ExpressionMatrix.",
           arg("directoryName"),
           arg("allowReadOnly")=false
       )

       // Get the total number of genes or cells currently in the system.
       .def("geneCount",
           &ExpressionMatrix::geneCount,
           "Return the total number of genes."
       )
       .def("cellCount",
           &ExpressionMatrix::cellCount,
           "Return the total number of cells."
       )

       // Genes.
       .def("addGene",
           &ExpressionMatrix::addGene,
           "Add a new gene with a given name.",
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
           (
               CellId (ExpressionMatrix::*)
               (const string& jsonString)
           )
           &ExpressionMatrix::addCell,
           "Add a cell described by a JSON string.",
           arg("jsonString")
       )
       .def("addCells",
           &ExpressionMatrix::addCells,
           "Add cells and their meta data from input files in delimited format.",
           arg("expressionCountsFileName"),
           arg("expressionCountsFileSeparators") = ",",
           arg("cellMetaDataFileName"),
           arg("cellMetaDataFileSeparators") = ","
       )
       .def("addCellsFromHdf5",
           &ExpressionMatrix::addCellsFromHdf5,
           "Add cells from an input file in hdf5 format.",
           arg("fileName"),
           arg("cellNamePrefix"),
           arg("cellMetaDataArgument"),
           arg("totalExpressionCountThreshold")
       )
       .def("addCellsFromBioHub1",
           &ExpressionMatrix::addCellsFromBioHub1,
           "Add cells in the format used by the BioHub pipeline, July 2017.",
           arg("expressionCountsFileName"),
           arg("initialMetaDataCount"),
           arg("finalMetaDataCount"),
           arg("plateMetaDataFileName")
       )
       .def("addCellsFromBioHub2",
           &ExpressionMatrix::addCellsFromBioHub2,
           "Add cells in the format used by the BioHub pipeline, September 2017.",
           arg("plateFileName"),
           arg("totalExpressionCountThreshold")
       )
       .def("addCellMetaData",
           &ExpressionMatrix::addCellMetaData,
           "Add meta data for existing cells from an input csv file.",
           arg("cellMetaDataFileName")
       )



       // Accessors for cell meta data.
       .def
       (
           "getCellMetaDataValue",
           (
               string (ExpressionMatrix::*)
               (CellId, const string& name) const
           )
           &ExpressionMatrix::getCellMetaData,
           "Get a single meta data value for a given cell id.",
           arg("cellId"),
           arg("name")
       )
       .def
       (
           "getCellMetaData",
           (
               vector< pair<string, string> > (ExpressionMatrix::*)
               (CellId) const
           )
           &ExpressionMatrix::getCellMetaData,
           "Get all meta data values for a given cell id.",
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
           "Get all meta data fields for a given list of cell ids.",
           arg("cellIds")
       )
       .def("cellIdFromString",
           &ExpressionMatrix::cellIdFromString,
           "Return the cellId corresponding to a given cell name (or cell id string).",
           arg("cellString")
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
           "Get the expression count for a given cell id and gene id.",
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
           "Get the expression count for a given cell id and gene name.",
           arg("cellId"),
           arg("geneName")
       )
       .def("getCellExpressionCounts",
           &ExpressionMatrix::getCellExpressionCounts,
           "Get all the expression counts for a given cell id.",
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
           "Get the expression counts for a given gene id and for a specified list of cell ids.",
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
           "Get the expression counts for a given gene name and for a specified list of cell ids.",
           arg("cellIds"),
           arg("geneName")
       )
       .def("getCellsExpressionCounts",
           &ExpressionMatrix::getCellsExpressionCounts,
           "Get all the expression counts for a specified list of cell ids.",
           arg("cellIds")
           )
       .def("getCellsExpressionCountsForGenes",
           &ExpressionMatrix::getCellsExpressionCountsForGenes,
           "Get all the expression counts for a specified list of cell ids and gene ids.",
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
           "Create a gene set using an information content criterion.",
           arg("existingGeneSetName"),
           arg("cellSetName"),
           arg("normalizationMethod"),
           arg("geneInformationContentThreshold"),
           arg("newGeneSetName")
       )
       .def("createGeneSetIntersection", &ExpressionMatrix::createGeneSetIntersection,
           "Create a new gene set as the intersection of two or more existing gene sets.",
           arg("inputSetsNames"),
           arg("outputSetName")
       )
       .def("createGeneSetUnion", &ExpressionMatrix::createGeneSetUnion,
           "Create a new gene set as the union of two or more existing gene sets.",
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
           "Create a new gene set as the set difference of two existing gene sets.",
           arg("inputSetName0"),
           arg("inputSetName1"),
           arg("oututSetName")
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
           "Create a new cell set using cell meta data.",
           arg("cellSetName"),
           arg("metaDataFieldName"),
           arg("matchString"),
           arg("useRegex")
       )
       .def("createCellSetUsingNumericMetaDataGreaterThan",
           &ExpressionMatrix::createCellSetUsingNumericMetaDataGreaterThan,
           "Create a new cell set consisting of cells for which a specified meta data field "
           "is numeric and greater than a specified value.",
           arg("cellSetName"),
           arg("metaDataFieldName"),
           arg("lowerBound")
       )
       .def("createCellSetUsingNumericMetaDataLessThan",
           &ExpressionMatrix::createCellSetUsingNumericMetaDataLessThan,
           "Create a new cell set consisting of cells for which a specified meta data field "
           "is numeric and less than a specified value.",
           arg("cellSetName"),
           arg("metaDataFieldName"),
           arg("upperBound")
       )
       .def("createCellSetUsingNumericMetaDataBetween",
           &ExpressionMatrix::createCellSetUsingNumericMetaDataBetween,
           "Create a new cell set consisting of cells for which a specified meta data field "
           "is numeric and in a specified interval.",
           arg("cellSetName"),
           arg("metaDataFieldName"),
           arg("lowerBound"),
           arg("upperBound")
       )
       .def("createCellSetIntersection",
           &ExpressionMatrix::createCellSetIntersection,
           "Create a new cell set as the intersection of two or more existing cell sets.",
           arg("inputSetsNames"),
           arg("outputSetName")
       )
       .def("createCellSetUnion",
           &ExpressionMatrix::createCellSetUnion,
           "Create a new cell set as the union of two or more existing cell sets.",
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
           "Create a new cell set as the set difference of two existing cell sets.",
           arg("inputSetName0"),
           arg("inputSetName1"),
           arg("outputSetName")
       )
       .def("getCellSetNames",
           &ExpressionMatrix::getCellSetNames,
           "Get the names of all existing cell sets."
       )
       .def("getCellSet",
           &ExpressionMatrix::getCellSet,
           "Get the ids of the cells in an existing cell set.",
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
           "Remove an existing cell set.",
           arg("cellSetName")
       )



       // Compute cell similarity.
       .def("computeCellSimilarity",
           &ExpressionMatrix::computeCellSimilarity,
           "Compute the similarity between two cells, taking into account all genes.",
           arg("cellId0"),
           arg("cellId1")
       )
       .def("computeApproximateLshCellSimilarity",
           &ExpressionMatrix::computeApproximateLshCellSimilarity,
           "For debugging/testing use only."
       )
       .def("writeLshSimilarityComparison",
           &ExpressionMatrix::writeLshSimilarityComparison,
           "For debugging/testing use only."
       )
       .def("writeLshSimilarityComparisonSlow",
           &ExpressionMatrix::writeLshSimilarityComparisonSlow,
           "For debugging/testing use only."
       )

       // Find similar pairs of cells.
       .def("findSimilarPairs0",
           &ExpressionMatrix::findSimilarPairs0,
           "Find pairs of similar cells using direct computation.",
           arg("geneSetName"),
           arg("cellSetName"),
           arg("similarPairsName"),
           arg("k"),
           arg("similarityThreshold")
       )
       .def("findSimilarPairs1",
           &ExpressionMatrix::findSimilarPairs1,
           "Obsolete. Use findSimilarPairs3 instead."
       )
       .def("findSimilarPairs2",
           &ExpressionMatrix::findSimilarPairs2,
           "Prototype code, not ready for production use."
       )
       .def("findSimilarPairs3",
           &ExpressionMatrix::findSimilarPairs3,
           "Find pairs of similar cells using Locality Sensitive Hashing, but still looping over all possible cell pairs.",
           arg("geneSetName") = "AllGenes",
           arg("cellSetName") = "AllCells",
           arg("similarPairsName"),
           arg("k") = 20,
           arg("similarityThreshold") = 0.5,
           arg("lshCount") = 1024,
           arg("seed") = 231
       )
       .def("writeSimilarPairs",
           &ExpressionMatrix::writeSimilarPairs,
           "For debugging/testing use only."
       )
       .def("analyzeSimilarPairs",
           &ExpressionMatrix::analyzeSimilarPairs,
           "For debugging/testing use only."
       )
       .def("analyzeLsh",
           &ExpressionMatrix::analyzeLsh,
           "For debugging/testing use only."
       )



       // Cell graphs.
       .def("getCellGraphNames",
           &ExpressionMatrix::getCellGraphNames,
           "Get the names of all currently defined cell graphs."
       )
       .def
       (
           "createCellGraph",
           (
               void (ExpressionMatrix::*)
               (const string&, const string&, const string&, double, size_t)
           )
           &ExpressionMatrix::createCellGraph,
           "Create a new cell graph.",
           arg("graphName"),
           arg("cellSetName") = "AllCells",
           arg("similarPairsName"),
           arg("similarityThreshold"),
           arg("k")
       )
       .def("computeCellGraphLayout",
           &ExpressionMatrix::computeCellGraphLayout,
           "Compute the two-dimensional layout of a cell graph to be used for visualization.",
           arg("graphName")
       )
       .def("getCellGraphVertices",
           &ExpressionMatrix::getCellGraphVertices,
           "Get the cell ids corresponding to the vertices of a cell graph.",
           arg("graphName")
       )
       .def("getCellGraphEdges",
           &ExpressionMatrix::getCellGraphEdges,
           "Get the pairs of cell ids corresponding to the edges of a cell graph.",
           arg("graphName")
       )



       // Cluster graphs.
       .def
       (
           "createClusterGraph",
           (
               void (ExpressionMatrix::*)
               (const string&, const ClusterGraphCreationParameters&, const string&)
           )
           &ExpressionMatrix::createClusterGraph,
           "Create a new cluster graph by running clustering on an existing cell graph.",
           arg("cellGraphName"),
           arg("clusterGraphCreationParameters"),
           arg("clusterGraphName")
       )
       .def("getClusterGraphVertices",
           &ExpressionMatrix::getClusterGraphVertices,
           "Get the cluster ids for the vertices of a cluster graph.",
           arg("clusterGraphName")
       )
       .def("getClusterGraphGenes",
           &ExpressionMatrix::getClusterGraphGenes,
           "Get the gene ids of the genes used by a cluster graph.",
           arg("clusterGraphName")
       )
       .def("getClusterCells",
           &ExpressionMatrix::getClusterCells,
           "Get the cell ids of a cluster.",
           arg("clusterGraphName"),
           arg("clusterId")
       )
       .def("getClusterAverageExpression",
           &ExpressionMatrix::getClusterAverageExpression,
           "Get L2-normalized average gene expression for a cluster.",
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
           "Create meta data from the cluster ids stored in a ClusterGraph.",
           arg("clusterGraphName"),
           arg("metaDataName")
       )



       // Run the http server.
       .def("explore",
           &ExpressionMatrix::explore,
           "Run the http server.",
           arg("serverParameters")
       )



       // Functions used only for testing or debugging.
       .def("testExpressionMatrixSubset",
           &ExpressionMatrix::testExpressionMatrixSubset,
           "For debugging/testing use only."
       )
       ;



    // Class ExpressionMatrixCreationParameters.
    class_<ExpressionMatrixCreationParameters>(
        module,
        "ExpressionMatrixCreationParameters",
        "Class used to store creation parameters for a new expression matrix.")
        .def(init<>())
        .def_readwrite("geneCapacity", &ExpressionMatrixCreationParameters::geneCapacity)
        .def_readwrite("cellCapacity", &ExpressionMatrixCreationParameters::cellCapacity)
        .def_readwrite("cellMetaDataNameCapacity", &ExpressionMatrixCreationParameters::cellMetaDataNameCapacity)
        .def_readwrite("cellMetaDataValueCapacity", &ExpressionMatrixCreationParameters::cellMetaDataValueCapacity)
        ;



    // Class ClusterGraphCreationParameters.
    class_<ClusterGraphCreationParameters>(module, "ClusterGraphCreationParameters")
        .def(init<>())
        .def_readwrite("stableIterationCount", &ClusterGraphCreationParameters::stableIterationCount)
        .def_readwrite("maxIterationCount", &ClusterGraphCreationParameters::maxIterationCount)
        .def_readwrite("seed", &ClusterGraphCreationParameters::seed)
        .def_readwrite("minClusterSize", &ClusterGraphCreationParameters::minClusterSize)
        .def_readwrite("maxConnectivity", &ClusterGraphCreationParameters::maxConnectivity)
        .def_readwrite("similarityThreshold", &ClusterGraphCreationParameters::similarityThreshold)
        ;



    // Class ServerParameters.
    class_<ServerParameters>(module, "ServerParameters")
        .def(init<>())
        .def_readwrite("port", &ServerParameters::port)
        .def_readwrite("docDirectory", &ServerParameters::docDirectory)
        ;



    // Class CellGraphVertexInfo.
    class_<CellGraphVertexInfo>(
        module,
        "CellGraphVertexInfo",
        "Class used to store information about a vertex of a cell graph.")
        .def(init<>())
        .def_readonly("cellId", &CellGraphVertexInfo::cellId)
        .def("x", &CellGraphVertexInfo::x)
        .def("y", &CellGraphVertexInfo::y)
        ;


    // Enum class NormalizationMethod.
    enum_<NormalizationMethod>(
        module,
        "NormalizationMethod",
        "Various ways to normalize gene expressions.")
        .value(normalizationMethodToShortString(NormalizationMethod::None).c_str(),       NormalizationMethod::None)
        .value(normalizationMethodToShortString(NormalizationMethod::L1).c_str(),         NormalizationMethod::L1)
        .value(normalizationMethodToShortString(NormalizationMethod::L2).c_str(),         NormalizationMethod::L2)
        .value(normalizationMethodToShortString(NormalizationMethod::Invalid).c_str(),    NormalizationMethod::Invalid)
        .export_values()
        ;


    // Some non-member functions used only for testing or debugging.
    module.def("testMemoryMappedVector",
        testMemoryMappedVector,
        "For debugging/testing use only."
        );
    module.def("testMemoryMappedVectorOfLists",
        testMemoryMappedVectorOfLists,
        "For debugging/testing use only."
        );
    module.def("testMemoryMappedStringTable",
        testMemoryMappedStringTable,
        "For debugging/testing use only."
        );
}
