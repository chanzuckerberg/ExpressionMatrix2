
// CZI.
#include "ExpressionMatrix.hpp"
#include "MemoryMappedVectorOfLists.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

// Boost libraries.
#include <boost/python.hpp>
using namespace boost::python;



BOOST_PYTHON_MODULE(ExpressionMatrix2)
{
    // Overloaded functions need special handling.
    // See http://www.boost.org/doc/libs/1_58_0/libs/python/doc/tutorial/doc/html/python/functions.html#python.overloading
    CellId (ExpressionMatrix::*addCell)(const string&, size_t)
        = &ExpressionMatrix::addCell;
    bool (ExpressionMatrix::*createCellSetUsingMetaData)(const string&, const string&, const string&)
        = &ExpressionMatrix::createCellSetUsingMetaData;
    string (ExpressionMatrix::*getMetaData)(CellId, const string& name) const
        = &ExpressionMatrix::getMetaData;



    // Expose class ExpressionMatrix to Python.
    class_<ExpressionMatrix, boost::noncopyable>("ExpressionMatrix", init<string, ExpressionMatrixCreationParameters>())
       .def(init<string>())
       .def("geneCount", &ExpressionMatrix::geneCount)
       .def("cellCount", &ExpressionMatrix::cellCount)
       .def("addGene", &ExpressionMatrix::addGene)
       .def("addCell", addCell)
       .def("addCells", &ExpressionMatrix::addCells)
       .def("addCellsFromHdf5", &ExpressionMatrix::addCellsFromHdf5)
       .def("addCellsFromBioHub", &ExpressionMatrix::addCellsFromBioHub)
       .def("addCellMetaData", &ExpressionMatrix::addCellMetaData)
       .def("getMetaData", getMetaData)
       .def("createCellSetUsingMetaData", createCellSetUsingMetaData)
       .def("createCellSetIntersection", &ExpressionMatrix::createCellSetIntersection)
       .def("createCellSetUnion", &ExpressionMatrix::createCellSetUnion)
       .def("computeCellSimilarity", &ExpressionMatrix::computeCellSimilarity)
       .def("computeApproximateCellSimilarity", &ExpressionMatrix::computeApproximateCellSimilarity)
       .def("computeApproximateLshCellSimilarity", &ExpressionMatrix::computeApproximateLshCellSimilarity)
       .def("writeLshSimilarityComparison", &ExpressionMatrix::writeLshSimilarityComparison)
       .def("writeLshSimilarityComparisonSlow", &ExpressionMatrix::writeLshSimilarityComparisonSlow)
       .def("analyzeAllPairs", &ExpressionMatrix::analyzeAllPairs)
       .def("findSimilarPairs0", &ExpressionMatrix::findSimilarPairs0)
       .def("findSimilarPairs1", &ExpressionMatrix::findSimilarPairs1)
       .def("findSimilarPairs2", &ExpressionMatrix::findSimilarPairs2)
       .def("writeSimilarPairs", &ExpressionMatrix::writeSimilarPairs)
       .def("createCellSimilarityGraph", &ExpressionMatrix::createCellSimilarityGraph)
       .def("explore", &ExpressionMatrix::explore)
       ;

    // Expose class ExpressionMatrixCreationParameters to Python.
    class_<ExpressionMatrixCreationParameters>("ExpressionMatrixCreationParameters", init<>())
        .def_readwrite("geneCapacity", &ExpressionMatrixCreationParameters::geneCapacity)
        .def_readwrite("cellCapacity", &ExpressionMatrixCreationParameters::cellCapacity)
        .def_readwrite("cellMetaDataNameCapacity", &ExpressionMatrixCreationParameters::cellMetaDataNameCapacity)
        .def_readwrite("cellMetaDataValueCapacity", &ExpressionMatrixCreationParameters::cellMetaDataValueCapacity)
        ;

    // Expose class ServerParameters to Python.
    class_<ServerParameters>("ServerParameters", init<>())
        .def_readwrite("port", &ServerParameters::port)
        .def_readwrite("docDirectory", &ServerParameters::docDirectory)
        ;

#if 0
    // Expose class ExpressionMatrix::ApproximateCellSimilarity to Python.
    // An object of this type is returned by function ExpressionMatrix::computeApproximateCellSimilarity.
    class_<ApproximateCellSimilarity>("ApproximateCellSimilarity", no_init)
        .def_readwrite("estimate", &ApproximateCellSimilarity::estimate)
        .def_readwrite("lowerBound1", &ApproximateCellSimilarity::lowerBound1)
        .def_readwrite("upperBound1", &ApproximateCellSimilarity::upperBound1)
        .def_readwrite("lowerBound2", &ApproximateCellSimilarity::lowerBound2)
        .def_readwrite("upperBound2", &ApproximateCellSimilarity::upperBound2)
        .def_readwrite("upperBound3", &ApproximateCellSimilarity::upperBound3)
        ;
#endif

    // Expose some other functions to Python.
    def("testMemoryMappedVectorOfLists", testMemoryMappedVectorOfLists);
    def("testMemoryMappedStringTable", testMemoryMappedStringTable);
}
