#include "ExpressionMatrix.hpp"
#include "color.hpp"
#include "GeneGraph.hpp"
#include "SimilarGenePairs.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/graph/iteration_macros.hpp>



// Check that a gene graph does not exist,
// and throw an exception if it does.
void ExpressionMatrix::checkGeneGraphDoesNotExist(
    const string& geneGraphName) const
{
    if(geneGraphs.find(geneGraphName) != geneGraphs.end()) {
        throw runtime_error("Gene graph " + geneGraphName + " already exists.");
    }

}

// Return a reference to the gene graph with a given name,
// and throw and exception if not found.
GeneGraph& ExpressionMatrix::getGeneGraph(const string& geneGraphName)
{
    const auto it = geneGraphs.find(geneGraphName);
    if(it == geneGraphs.end()) {
        throw runtime_error("Gene graph " + geneGraphName + " does not exists.");
    }
    return *(it->second);
}



void ExpressionMatrix::createGeneGraph(
    const string& geneGraphName,
    const string& geneSetName,
    const string& similarGenePairsName,
    int k,
    double similarityThreshold)
{
    createGeneGraph(cout,
        geneGraphName, geneSetName, similarGenePairsName, k, similarityThreshold);
}
void ExpressionMatrix::createGeneGraph(
    ostream& out,
    const string& geneGraphName,
    const string& geneSetName,
    const string& similarGenePairsName,
    int k,
    double similarityThreshold)
{
    // Check that a signature graph with this name does not already exist.
    checkSignatureGraphDoesNotExist(geneGraphName);

    // Locate the gene set and verify that it is not empty.
    const auto& it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    GeneSet& geneSet = it->second;
    const GeneId geneCount = GeneId(geneSet.size());
    if(geneCount == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Access the similar gene pairs.
    SimilarGenePairs similarGenePairs(directoryName, similarGenePairsName , true);

    // Now we have everything we need. Create the gene graph.
    const shared_ptr<GeneGraph> geneGraphPointer =
        make_shared<GeneGraph>(
        out,
        geneSet,
        directoryName,
        similarGenePairsName,
        similarityThreshold,
        k);
    geneGraphs.insert(make_pair(geneGraphName, geneGraphPointer));
}



void ExpressionMatrix::removeGeneGraph(const string& geneGraphName)
{
    const auto it = geneGraphs.find(geneGraphName);
    if(it == geneGraphs.end()) {
        throw runtime_error("Gene graph " + geneGraphName + " does not exists.");
    }
    geneGraphs.erase(it);
}



void ExpressionMatrix::colorGeneGraphBlack(GeneGraph& geneGraph) const
{
    BGL_FORALL_VERTICES(v, geneGraph, GeneGraph) {
        geneGraph[v].color= "black";
    }
}



void ExpressionMatrix::colorGeneGraphByMetaData(GeneGraph& geneGraph, const string& metaDataName) const
{

    // Count the number of vertices with each value of the given meta data field.
    map<string, GeneId> frequency;
    BGL_FORALL_VERTICES(v, geneGraph, GeneGraph) {
        const GeneGraphVertex& vertex = geneGraph[v];
        const GeneId globalGeneId = vertex.globalGeneId;
        const string metaDataValue = getGeneMetaData(globalGeneId, metaDataName);
        const auto it = frequency.find(metaDataValue);
        if(it == frequency.end()) {
            frequency.insert(make_pair(metaDataValue, 1));
        } else {
            ++(it->second);
        }
    }

    // Sort by decreasing frequency, without counting empty values.
    vector< pair<string, GeneId> > sortedFrequency;
    for(const auto& p: frequency) {
        if(!p.first.empty()) {
            sortedFrequency.push_back(p);
        }
    }
    sort(sortedFrequency.begin(), sortedFrequency.end(), OrderPairsBySecondGreater< pair<string, GeneId> >());



    // Assign colors to meta data values.
    // Vertices for which the value is empty are colored black.
    // The remaining values are assigned colors from colorPalette1
    // in order of decreasing frequency.
    map<string, string> colorMap;
    for(size_t i=0; (i<sortedFrequency.size()) && (i<12); i++) {
        const auto& p = sortedFrequency[i];
        colorMap[p.first] = colorPalette1(i);
    }

    // Assign colors to the vertices.
    BGL_FORALL_VERTICES(v, geneGraph, GeneGraph) {
        const GeneGraphVertex& vertex = geneGraph[v];
        const GeneId globalGeneId = vertex.globalGeneId;
        const string metaDataValue = getGeneMetaData(globalGeneId, metaDataName);
        const auto it = colorMap.find(metaDataValue);
        if(it == colorMap.end()) {
            geneGraph[v].color= "black";
        } else {
            geneGraph[v].color = it->second;
        }
    }
}

