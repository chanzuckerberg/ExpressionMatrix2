#include "SimilarGenePairs.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



// Create a new SimilarGenePairs object.
SimilarGenePairs::SimilarGenePairs(
    const string& name,
    size_t k,
    const GeneSet& geneSetArgument,
    const vector< vector<Pair> >& pairsArgument)
{
    const GeneId geneCount = GeneId(geneSetArgument.size());

    info.createNew(name + "-Info");
    info->k = k;
    info->geneCount = geneCount;

    pairs.createNew(name + "-Pairs", k*size_t(geneCount));

    // Initialize the geneInfo vector.
    geneInfo.createNew(name + "-GeneInfo", geneCount);

    // Make a copy of the gene set. The copy is owned by the SimilarGenePairs object.
    geneSetArgument.makeCopy(geneSet, name + "-GeneSet");

    // Store the pairs.
    for(GeneId geneId=0; geneId!=geneCount; geneId++) {
        const auto& v = pairsArgument[geneId];
        const GeneId n = GeneId(v.size());
        CZI_ASSERT(n <= k);
        geneInfo[geneId].usedCount = n;
        copy(v.begin(), v.end(), begin(geneId));
    }

}

// Access an existing SimilarPairs object.
SimilarGenePairs::SimilarGenePairs(const string& name, bool allowReadOnly)
{
    info.accessExistingReadOnly(name + "-Info");
    geneInfo.accessExistingReadOnly(name + "-GeneInfo");
    pairs.accessExistingReadOnly(name + "-Pairs");
    geneSet.accessExisting(name + "-GeneSet", allowReadOnly);
}



void SimilarGenePairs::remove()
{
    info.remove();
    geneInfo.remove();
    pairs.remove();
    geneSet.remove();
}
