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

    // Class ExpressionMatrix.
    class_<ExpressionMatrix>(module, "ExpressionMatrix")
       .def(init<string, ExpressionMatrixCreationParameters>())
       .def(init<string, bool>())

       // Get the total number of genes or cells currently in the system.
       .def("geneCount", &ExpressionMatrix::geneCount)
       .def("cellCount", &ExpressionMatrix::cellCount)

       // Genes
       .def("addGene", &ExpressionMatrix::addGene)
       .def("geneName", &ExpressionMatrix::geneName)
       .def("geneIdFromName", &ExpressionMatrix::geneIdFromName)



       // Various ways to add cells.
       .def
       (
           "addCell",
           (
               CellId (ExpressionMatrix::*)
               (const string& jsonString)
           )
           &ExpressionMatrix::addCell
       )
       .def("addCells", &ExpressionMatrix::addCells)
       .def("addCellsFromHdf5", &ExpressionMatrix::addCellsFromHdf5)
       .def("addCellsFromBioHub1", &ExpressionMatrix::addCellsFromBioHub1)
       .def("addCellsFromBioHub2", &ExpressionMatrix::addCellsFromBioHub2)
       .def("addCellMetaData", &ExpressionMatrix::addCellMetaData)



       // Accessors for cell meta data.
       .def
       (
           "getCellMetaDataValue",
           (
               string (ExpressionMatrix::*)
               (CellId, const string& name) const
           )
           &ExpressionMatrix::getCellMetaData
       )
       .def
       (
           "getCellMetaData",
           (
               vector< pair<string, string> > (ExpressionMatrix::*)
               (CellId) const
           )
           &ExpressionMatrix::getCellMetaData
       )
       .def
       (
           "getCellsMetaData",
           (
               vector< vector< pair<string, string> > > (ExpressionMatrix::*)
               (const vector<CellId>&) const
           )
           &ExpressionMatrix::getCellMetaData
       )
       .def("cellIdFromString", &ExpressionMatrix::cellIdFromString)



       // Accessors for expression counts.
       .def
       (
           "getCellExpressionCount",
           (
               float (ExpressionMatrix::*)
               (CellId, GeneId) const
           )
           &ExpressionMatrix::getCellExpressionCount
       )
       .def
       (
           "getCellExpressionCountFromGeneName",
           (
               float (ExpressionMatrix::*)
               (CellId, const string& geneName) const
           )
           &ExpressionMatrix::getCellExpressionCount
       )
       .def("getCellExpressionCounts", &ExpressionMatrix::getCellExpressionCounts)
       .def
       (
           "getCellsExpressionCount",
           (
               vector<float> (ExpressionMatrix::*)
               (const vector<CellId>&, GeneId) const
           )
           &ExpressionMatrix::getCellsExpressionCount
       )
       .def
       (
           "getCellsExpressionCountFromGeneName",
           (
               vector<float> (ExpressionMatrix::*)
               (const vector<CellId>&, const string& geneName) const
           )
           &ExpressionMatrix::getCellsExpressionCount
       )
       .def("getCellsExpressionCounts", &ExpressionMatrix::getCellsExpressionCounts)
       .def("getCellsExpressionCountsForGenes", &ExpressionMatrix::getCellsExpressionCountsForGenes)



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
           &ExpressionMatrix::createGeneSetUsingInformationContent
       )
       .def("createGeneSetIntersection", &ExpressionMatrix::createGeneSetIntersection)
       .def("createGeneSetUnion", &ExpressionMatrix::createGeneSetUnion)
       .def
       (
           "createGeneSetDifference",
           (
               bool (ExpressionMatrix::*)
               (const string&, const string&, const string&)
           )
           &ExpressionMatrix::createGeneSetDifference
       )



        // Cell sets.
       .def
       (
           "createCellSetUsingMetaData",
           (
               bool (ExpressionMatrix::*)
               (const string&, const string&, const string&, bool)
           )
           &ExpressionMatrix::createCellSetUsingMetaData
       )
       .def("createCellSetUsingNumericMetaDataGreaterThan", &ExpressionMatrix::createCellSetUsingNumericMetaDataGreaterThan)
       .def("createCellSetUsingNumericMetaDataLessThan", &ExpressionMatrix::createCellSetUsingNumericMetaDataLessThan)
       .def("createCellSetUsingNumericMetaDataBetween", &ExpressionMatrix::createCellSetUsingNumericMetaDataBetween)
       .def("createCellSetIntersection", &ExpressionMatrix::createCellSetIntersection)
       .def("createCellSetUnion", &ExpressionMatrix::createCellSetUnion)
       .def
       (
           "createCellSetDifference",
           (
               bool (ExpressionMatrix::*)
               (const string&, const string&, const string&)
           )
           &ExpressionMatrix::createCellSetDifference
       )
       .def("getCellSetNames", &ExpressionMatrix::getCellSetNames)
       .def("getCellSet", &ExpressionMatrix::getCellSet)
       .def
       (
           "removeCellSet",
           (
               void (ExpressionMatrix::*)
               (const string&)
           )
           &ExpressionMatrix::removeCellSet
       )



       // Compute cell similarity.
       .def("computeCellSimilarity", &ExpressionMatrix::computeCellSimilarity)
       .def("computeApproximateLshCellSimilarity", &ExpressionMatrix::computeApproximateLshCellSimilarity)
       .def("writeLshSimilarityComparison", &ExpressionMatrix::writeLshSimilarityComparison)
       .def("writeLshSimilarityComparisonSlow", &ExpressionMatrix::writeLshSimilarityComparisonSlow)

       // Find similar pairs of cells.
       .def("findSimilarPairs0", &ExpressionMatrix::findSimilarPairs0)
       .def("findSimilarPairs1", &ExpressionMatrix::findSimilarPairs1)
       .def("findSimilarPairs2", &ExpressionMatrix::findSimilarPairs2)
       .def("findSimilarPairs3", &ExpressionMatrix::findSimilarPairs3)
       .def("writeSimilarPairs", &ExpressionMatrix::writeSimilarPairs)
       .def("analyzeSimilarPairs", &ExpressionMatrix::analyzeSimilarPairs)
       .def("analyzeLsh", &ExpressionMatrix::analyzeLsh)

       // Cell graphs.
       .def("getCellGraphNames", &ExpressionMatrix::getCellGraphNames)
       .def
       (
           "createCellGraph",
           (
               void (ExpressionMatrix::*)
               (const string&, const string&, const string&, double, size_t)
           )
           &ExpressionMatrix::createCellGraph
       )
       .def("computeCellGraphLayout", &ExpressionMatrix::computeCellGraphLayout)
       .def("getCellGraphVertices", &ExpressionMatrix::getCellGraphVertices)
       .def("getCellGraphEdges", &ExpressionMatrix::getCellGraphEdges)

       // Cluster graphs.
       .def
       (
           "createClusterGraph",
           (
               void (ExpressionMatrix::*)
               (const string&, const ClusterGraphCreationParameters&, const string&)
           )
           &ExpressionMatrix::createClusterGraph
       )
       .def("getClusterGraphVertices", &ExpressionMatrix::getClusterGraphVertices)
       .def("getClusterGraphGenes", &ExpressionMatrix::getClusterGraphGenes)
       .def("getClusterCells", &ExpressionMatrix::getClusterCells)
       .def("getClusterAverageExpression", &ExpressionMatrix::getClusterAverageExpression)
       .def
       (
           "createMetaDataFromClusterGraph",
           (
               void (ExpressionMatrix::*)
               (const string&, const string&)
           )
           &ExpressionMatrix::createMetaDataFromClusterGraph
       )

       // Run the http server.
       .def("explore", &ExpressionMatrix::explore)

       // Functions used only for testing or debugging.
       .def("testExpressionMatrixSubset", &ExpressionMatrix::testExpressionMatrixSubset)
       ;



    // Class ExpressionMatrixCreationParameters.
    class_<ExpressionMatrixCreationParameters>(module,"ExpressionMatrixCreationParameters")
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
        ;



    // Class ServerParameters.
    class_<ServerParameters>(module, "ServerParameters")
        .def(init<>())
        .def_readwrite("port", &ServerParameters::port)
        .def_readwrite("docDirectory", &ServerParameters::docDirectory)
        ;



    // Class CellGraphVertexInfo.
    class_<CellGraphVertexInfo>(module, "CellGraphVertexInfo")
        .def(init<>())
        .def_readonly("cellId", &CellGraphVertexInfo::cellId)
        .def("x", &CellGraphVertexInfo::x)
        .def("y", &CellGraphVertexInfo::y)
        ;


    // Enum class NormalizationMethod.
    enum_<NormalizationMethod>(module, "NormalizationMethod")
      .value(normalizationMethodToShortString(NormalizationMethod::None).c_str(),       NormalizationMethod::None)
      .value(normalizationMethodToShortString(NormalizationMethod::L1).c_str(),         NormalizationMethod::L1)
      .value(normalizationMethodToShortString(NormalizationMethod::L2).c_str(),         NormalizationMethod::L2)
      .value(normalizationMethodToShortString(NormalizationMethod::Invalid).c_str(),    NormalizationMethod::Invalid)
      .export_values()
      ;


    // Some non-member functions used only for testing or debugging.
    module.def("testMemoryMappedVector", testMemoryMappedVector);
    module.def("testMemoryMappedVectorOfLists", testMemoryMappedVectorOfLists);
    module.def("testMemoryMappedStringTable", testMemoryMappedStringTable);
}
