#include "SimilarGenePairs.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



// Create a new SimilarGenePairs object.
SimilarGenePairs::SimilarGenePairs(
    const string& directoryName,
    const string& similarGenePairsName,
    const string& geneSetName,
    const string& cellSetName,
    size_t k,
    const vector< vector<Pair> >& pairsArgument)
{
    // Access the gene set and the cell set.
    accessGeneSet(directoryName, geneSetName);
    accessCellSet(directoryName, cellSetName);
    if(pairsArgument.size() != geneSet.size()) {
        throw runtime_error("The gene pairs vector to create similar gene pairs " +
            similarGenePairsName +
            " has length inconsistent with gene set " + geneSetName
            );
    }

    // Create the info object and fill it in.
    const string pathBaseName = getPathBaseName(directoryName, similarGenePairsName);
    info.createNew(pathBaseName + "-Info");
    info->k = k;
    info->geneSetName = geneSetName;
    info->geneSetHash = geneSet.genes().hash();
    info->cellSetName = cellSetName;
    info->cellSetHash = cellSet.hash();

    // Create the remaining objects.
    pairs.createNew(pathBaseName + "-Pairs", k*size_t(geneSet.size()));
    geneInfo.createNew(pathBaseName + "-GeneInfo", geneSet.size());

    // Store the pairs that were passed in as the argument.
    for(GeneId localGeneId=0; localGeneId<geneSet.size(); localGeneId++){
        const vector<Pair>& pairsThisGene = pairsArgument[localGeneId];
        geneInfo[localGeneId].usedCount = GeneId(pairsThisGene.size());
        copy(pairsThisGene.begin(), pairsThisGene.end(), begin(localGeneId));
     }

}



// Access an existing SimilarPairs object.
SimilarGenePairs::SimilarGenePairs(
    const string& directoryName,
    const string& similarGenePairsName,
    bool allowReadOnly)
{
    // Access the info object.
    const string pathBaseName = getPathBaseName(directoryName, similarGenePairsName);
    info.accessExistingReadOnly(pathBaseName + "-Info");

    // Access the gene set and the cell set.
    accessGeneSet(directoryName, info->geneSetName);
    accessCellSet(directoryName, info->cellSetName);

    // Access the remaining objects.
    pairs.accessExistingReadOnly(pathBaseName + "-Pairs");
    geneInfo.accessExistingReadOnly(pathBaseName + "-GeneInfo");

    // Check that all is good.
    if(geneSet.genes().hash() != info->geneSetHash) {
        throw runtime_error("Hash for gene set " + string(info->geneSetName) +
            " is not consistent with the value at the time SimilarGenePairs object " +
            similarGenePairsName + " was created.");
    }
    if(cellSet.hash() != info->cellSetHash) {
        throw runtime_error("Hash for cell set " + string(info->cellSetName) +
            " is not consistent with the value at the time SimilarGenePairs object " +
            similarGenePairsName + " was created.");
    }
    if(geneInfo.size() != geneSet.size()) {
        throw runtime_error("SimilarGenePairs object " +
            similarGenePairsName + " has geneInfo vector of inconsistent length.");
    }
    if(pairs.size() != info->k * size_t(geneSet.size())) {
        throw runtime_error("SimilarGenePairs object " +
            similarGenePairsName + " has pairs vector of inconsistent length.");
    }
}

string SimilarGenePairs::getPathBaseName(
    const string& directoryName,
    const string& similarGenePairsName
    )
{
    return directoryName + "/SimilarGenePairs-" + similarGenePairsName;
}



void SimilarGenePairs::accessGeneSet(
    const string& directoryName,
    const string& geneSetName)
{
    geneSet.accessExisting(directoryName + "/GeneSet-" + geneSetName, true);
    if(!std::is_sorted(geneSet.genes().begin(), geneSet.genes().end())) {
        throw runtime_error("Gene set " + geneSetName + " is not sorted.");
    }
}
void SimilarGenePairs::accessCellSet(
    const string& directoryName,
    const string& cellSetName)
{
    cellSet.accessExisting(directoryName + "/CellSet-" + cellSetName, true);
    if(!std::is_sorted(cellSet.begin(), cellSet.end())) {
        throw runtime_error("Cell set " + cellSetName + " is not sorted.");
    }
}


void SimilarGenePairs::remove()
{
    info.remove();
    geneInfo.remove();
    pairs.remove();
    geneSet.remove();
}
