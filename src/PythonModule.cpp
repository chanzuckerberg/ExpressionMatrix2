// This file uses the Boost Python library to export portions of
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

// Boost libraries.
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;



namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        // Function used to expose an std::pair to Python.
        template<class X, class Y> void exposePair(const string& name)
        {
            class_< pair<X, Y> >(name.c_str(), init<>())    // Allow default construction
                .def(init<X, Y>())                          // Allow construction from an X and Y
                .def_readwrite("first",  &pair<X, Y>::first)
                .def_readwrite("second", &pair<X, Y>::second);
        }

        // Function used to expose an std::vector to Python.
        template<class T> void exposeVector(const string& name)
        {
            class_< vector<T> >(name.c_str())
                .def(vector_indexing_suite< vector<T> >());
        }
    }
}



BOOST_PYTHON_MODULE(ExpressionMatrix2)
{

    // Pair classes. These are instantiations of std::pair exposed to Python.
    // The two items in the pair can be accessed in Python as p.first and p.second.
    exposePair<int, int>("IntIntPair");
    exposePair<uint32_t, uint32_t>("UintUintPair");
    exposePair<string, string>("StringStringPair");
    exposePair<uint32_t, float>("UintFloatPair");



    // Python aliases for pair classes.
    scope().attr("IntPair") = scope().attr("IntIntPair");
    scope().attr("UintPair") = scope().attr("UintUintPair");
    scope().attr("CellIdCellIdPair") = scope().attr("UintUintPair");
    scope().attr("CellIdPair") = scope().attr("UintUintPair");
    scope().attr("StringPair") = scope().attr("StringStringPair");
    scope().attr("NameValuePair") = scope().attr("StringStringPair");
    scope().attr("GeneIdFloatPair") = scope().attr("UintFloatPair");
    scope().attr("ExpressionCount") = scope().attr("UintFloatPair");



    // Container classes. These are standard C++ vectors
    // exposed to Python using the vector_indexing_suite capability.
    // They have an API very similar to that of a Python list,
    // which means, for example, that they can be iterated on using typical
    // Python patterns. For example, if v is one of these containers,
    // the following Python code works with the expected semantics:
    // Iteration:
    //     for x in v:
    // Element access:
    //     x = v[i]
    //     v[i] = x
    // Number of elements in the container:
    //     n = len(v)
    // Convert to a Python list.
    //     vList = [x for x in v]
    // See http://www.boost.org/doc/libs/1_58_0/libs/python/doc/v2/indexing.html
    // for more information.
    // Because of the similarity with Python lists, we use Python names
    // that end with "List" even though these types map to C++ std::vector.
    exposeVector<int>("IntList");
    exposeVector<uint32_t>("UintList"); // Used for vector<GeneId> and vector<CellId>
    exposeVector< pair<int, int> >("IntIntPairList");
    exposeVector< pair<uint32_t, uint32_t> >("UintUintPairList");
    exposeVector<string>("StringList");
    exposeVector< pair<string, string> >("StringStringPairList");
    exposeVector< vector< pair<string, string> > >("StringStringPairListList");
    exposeVector<float>("FloatList");
    exposeVector< pair<uint32_t, float> >("UintFloatPairList");
    exposeVector< vector< pair<uint32_t, float> > >("UintFloatPairListList");
    exposeVector<CellGraphVertexInfo>("CellGraphVertexInfoList");



    // Python aliases for container types.

    // Python aliases for UintList.
    scope().attr("GeneIdList") = scope().attr("UintList");
    scope().attr("CellIdList") = scope().attr("UintList");

    // Python aliases for UintUintPairList.
    scope().attr("UintPairList") = scope().attr("UintUintPairList");
    scope().attr("CellIdCellIdPairList") = scope().attr("UintUintPairList");
    scope().attr("CellIdPairList") = scope().attr("UintUintPairList");

    // Python aliases for StringStringPairList.
    scope().attr("StringPairList") = scope().attr("StringStringPairList");
    scope().attr("NameValuePairList") = scope().attr("StringStringPairList");
    scope().attr("CellMetaData") = scope().attr("StringStringPairList");

    // Python aliases for StringStringPairListList.
    scope().attr("StringPairListList") = scope().attr("StringStringPairListList");
    scope().attr("NameValuePairListList") = scope().attr("StringStringPairListList");
    scope().attr("CellMetaDataList") = scope().attr("StringStringPairListList");

    // Python aliases for UintFloatPairList.
    scope().attr("GeneIdFloatPairList") = scope().attr("UintFloatPairList");
    scope().attr("ExpressionCountList") = scope().attr("UintFloatPairList");

    // Python aliases for UintFloatPairListList.
    scope().attr("GeneIdFloatPairListList") = scope().attr("UintFloatPairListList");
    scope().attr("ExpressionCountListList") = scope().attr("UintFloatPairListList");





    // Constants for an invalid gene id or cell id.
    scope().attr("invalidGeneId") = invalidGeneId;
    scope().attr("invalidCellId") = invalidCellId;



    // In the definitions below, overloaded functions need special handling.
    // See http://www.boost.org/doc/libs/1_58_0/libs/python/doc/tutorial/doc/html/python/functions.html#python.overloading



    // Class ExpressionMatrix.
    class_<ExpressionMatrix, boost::noncopyable>("ExpressionMatrix", init<string, ExpressionMatrixCreationParameters>())
       .def(init<string>())

       // Get the total number of genes or cells currently in the system.
       .def("geneCount", &ExpressionMatrix::geneCount)
       .def("cellCount", &ExpressionMatrix::cellCount)

       // Genes
       .def("addGene", &ExpressionMatrix::addGene)
       .def("geneName", &ExpressionMatrix::geneName)



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
               (const string&, const string&, const string&)
           )
           &ExpressionMatrix::createCellSetDifference
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



       // Compute cell similarity.
       .def("computeCellSimilarity", &ExpressionMatrix::computeCellSimilarity)
       .def("computeApproximateLshCellSimilarity", &ExpressionMatrix::computeApproximateLshCellSimilarity)
       .def("writeLshSimilarityComparison", &ExpressionMatrix::writeLshSimilarityComparison)
       .def("writeLshSimilarityComparisonSlow", &ExpressionMatrix::writeLshSimilarityComparisonSlow)

       // Find similar pairs of cells.
       .def("findSimilarPairs0", &ExpressionMatrix::findSimilarPairs0)
       .def("findSimilarPairs1", &ExpressionMatrix::findSimilarPairs1)
       .def("findSimilarPairs2", &ExpressionMatrix::findSimilarPairs2)
       .def("writeSimilarPairs", &ExpressionMatrix::writeSimilarPairs)
       .def("analyzeSimilarPairs", &ExpressionMatrix::analyzeSimilarPairs)

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
       .def("computeClusterGraphLayout", &ExpressionMatrix::computeClusterGraphLayout)

       // Run the http server.
       .def("explore", &ExpressionMatrix::explore)

       // Functions used only for testing or debugging.
       .def("testExpressionMatrixSubset", &ExpressionMatrix::testExpressionMatrixSubset)
       ;



    // Class ExpressionMatrixCreationParameters.
    class_<ExpressionMatrixCreationParameters>("ExpressionMatrixCreationParameters", init<>())
        .def_readwrite("geneCapacity", &ExpressionMatrixCreationParameters::geneCapacity)
        .def_readwrite("cellCapacity", &ExpressionMatrixCreationParameters::cellCapacity)
        .def_readwrite("cellMetaDataNameCapacity", &ExpressionMatrixCreationParameters::cellMetaDataNameCapacity)
        .def_readwrite("cellMetaDataValueCapacity", &ExpressionMatrixCreationParameters::cellMetaDataValueCapacity)
        ;

    // Class ClusterGraphCreationParameters.
    class_<ClusterGraphCreationParameters>("ClusterGraphCreationParameters", init<>())
        .def_readwrite("stableIterationCount", &ClusterGraphCreationParameters::stableIterationCount)
        .def_readwrite("maxIterationCount", &ClusterGraphCreationParameters::maxIterationCount)
        .def_readwrite("seed", &ClusterGraphCreationParameters::seed)
        ;



    // Class ServerParameters.
    class_<ServerParameters>("ServerParameters", init<>())
        .def_readwrite("port", &ServerParameters::port)
        .def_readwrite("docDirectory", &ServerParameters::docDirectory)
        ;

    // Class CellGraphVertexInfo.
    class_<CellGraphVertexInfo>("CellGraphVertexInfo", init<>())
        .def_readonly("cellId", &CellGraphVertexInfo::cellId)
        .def("x", &CellGraphVertexInfo::x)
        .def("y", &CellGraphVertexInfo::y)
        ;


    // Enum class NormalizationMethod.
    enum_<NormalizationMethod>("NormalizationMethod")
      .value(normalizationMethodToShortString(NormalizationMethod::None).c_str(),       NormalizationMethod::None)
      .value(normalizationMethodToShortString(NormalizationMethod::L1).c_str(),         NormalizationMethod::L1)
      .value(normalizationMethodToShortString(NormalizationMethod::L2).c_str(),         NormalizationMethod::L2)
      .value(normalizationMethodToShortString(NormalizationMethod::Invalid).c_str(),    NormalizationMethod::Invalid)
      .export_values()
      ;

    // Some non-member functions used only for testing or debugging.
    def("testMemoryMappedVector", testMemoryMappedVector);
    def("testMemoryMappedVectorOfLists", testMemoryMappedVectorOfLists);
    def("testMemoryMappedStringTable", testMemoryMappedStringTable);
}
