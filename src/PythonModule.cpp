// This file uses the Boost Python library to export portions of
// the C++ API to Python.
// If something is added or changed here, corresponding documentation
// changes should also be made in doc/PythonApiReference.html.


// CZI.
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
    exposePair<string, string>("StringStringPair");



    // Python aliases for pair classes.
    scope().attr("IntPair") = scope().attr("IntIntPair");
    scope().attr("StringPair") = scope().attr("StringStringPair");
    scope().attr("NameValuePair") = scope().attr("StringStringPair");



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
    exposeVector<string>("StringList");
    exposeVector< pair<string, string> >("StringStringPairList");
    exposeVector< vector< pair<string, string> > >("StringStringPairListList");
    exposeVector<CellGraphVertexInfo>("CellGraphVertexInfoList");



    // Python aliases for container types.

    // Python aliases for UintList.
    scope().attr("GeneIdList") = scope().attr("UintList");
    scope().attr("CellIdList") = scope().attr("UintList");

    // Python aliases for StringStringPairList.
    scope().attr("StringPairList") = scope().attr("StringStringPairList");
    scope().attr("NameValuePairList") = scope().attr("StringStringPairList");
    scope().attr("CellMetaData") = scope().attr("StringStringPairList");

    // Python aliases for StringStringPairListList.
    scope().attr("StringPairListList") = scope().attr("StringStringPairListList");
    scope().attr("NameValuePairListList") = scope().attr("StringStringPairListList");
    scope().attr("CellMetaDataList") = scope().attr("StringStringPairListList");




    // Constants for an invalid gene id or cell id.
    scope().attr("invalidGeneId") = invalidGeneId;
    scope().attr("invalidCellId") = invalidCellId;


    // Overloaded functions need special handling.
    // See http://www.boost.org/doc/libs/1_58_0/libs/python/doc/tutorial/doc/html/python/functions.html#python.overloading
    CellId (ExpressionMatrix::*addCell)(const string&)
        = &ExpressionMatrix::addCell;
    bool (ExpressionMatrix::*createCellSetUsingMetaData)(const string&, const string&, const string&)
        = &ExpressionMatrix::createCellSetUsingMetaData;
    string (ExpressionMatrix::*getCellMetaDataValue)(CellId, const string& name) const
        = &ExpressionMatrix::getCellMetaData;
    vector< pair<string, string> > (ExpressionMatrix::*getCellMetaData)(CellId) const
        = &ExpressionMatrix::getCellMetaData;
    vector< vector< pair<string, string> > > (ExpressionMatrix::*getCellsMetaData)(const vector<CellId>&) const
        = &ExpressionMatrix::getCellMetaData;
    void (ExpressionMatrix::*createGeneSetUsingInformationContent) (
        const string& existingGeneSetName,
        const string& cellSetName,
        NormalizationMethod normalizationMethod,
        double geneInformationContentThreshold,
        const string& newGeneSetName)
        = &ExpressionMatrix::createGeneSetUsingInformationContent;
    bool (ExpressionMatrix::*createGeneSetDifference) (
        const string& inputSet0,
        const string& inputSet1,
        const string& outputSet0)
        = &ExpressionMatrix::createGeneSetDifference;
    bool (ExpressionMatrix::*createCellSetDifference) (
        const string& inputSet0,
        const string& inputSet1,
        const string& outputSet0)
        = &ExpressionMatrix::createCellSetDifference;



    // Class ExpressionMatrix.
    class_<ExpressionMatrix, boost::noncopyable>("ExpressionMatrix", init<string, ExpressionMatrixCreationParameters>())
       .def(init<string>())

       // Get the total number of genes or cells currently in the system.
       .def("geneCount", &ExpressionMatrix::geneCount)
       .def("cellCount", &ExpressionMatrix::cellCount)

       // Add a gene.
       .def("addGene", &ExpressionMatrix::addGene)

       // Various ways to add cells.
       .def("addCell", addCell)
       .def("addCells", &ExpressionMatrix::addCells)
       .def("addCellsFromHdf5", &ExpressionMatrix::addCellsFromHdf5)
       .def("addCellsFromBioHub", &ExpressionMatrix::addCellsFromBioHub)
       .def("addCellMetaData", &ExpressionMatrix::addCellMetaData)

       // Get cell meta data.
       .def("getCellMetaDataValue", getCellMetaDataValue)
       .def("getCellMetaData", getCellMetaData)
       .def("getCellsMetaData", getCellsMetaData)
       .def("cellIdFromString", &ExpressionMatrix::cellIdFromString)

       // Gene sets.
       .def("createGeneSetUsingInformationContent", createGeneSetUsingInformationContent)
       .def("createGeneSetIntersection", &ExpressionMatrix::createGeneSetIntersection)
       .def("createGeneSetUnion", &ExpressionMatrix::createGeneSetUnion)
       .def("createGeneSetDifference", createGeneSetDifference)

       // Cell sets.
       .def("createCellSetUsingMetaData", createCellSetUsingMetaData)
       .def("createCellSetIntersection", &ExpressionMatrix::createCellSetIntersection)
       .def("createCellSetUnion", &ExpressionMatrix::createCellSetUnion)
       .def("createCellSetDifference", createCellSetDifference)
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
       .def("createCellGraph", &ExpressionMatrix::createCellGraph)
       .def("computeCellGraphLayout", &ExpressionMatrix::computeCellGraphLayout)
       .def("getCellGraphVertices", &ExpressionMatrix::getCellGraphVertices)

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
