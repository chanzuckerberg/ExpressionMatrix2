#include "ExpressionMatrix.hpp"
#include "GeneGraph.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



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
    checkSignatureGraphDoesNotExist(geneGraphName);

    // Locate the gene set and verify that it is not empty.
    const auto& it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = it->second;
    const GeneId geneCount = GeneId(geneSet.size());
    if(geneCount == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    CZI_ASSERT(0);
}



void ExpressionMatrix::removeGeneGraph(const string& geneGraphName)
{
    const auto it = geneGraphs.find(geneGraphName);
    if(it == geneGraphs.end()) {
        throw runtime_error("Gene graph " + geneGraphName + " does not exists.");
    }
    geneGraphs.erase(it);
}

