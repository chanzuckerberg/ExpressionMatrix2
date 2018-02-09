#include "ExpressionMatrix.hpp"
#include "ExpressionMatrixSubset.hpp"
#include "heap.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



void ExpressionMatrix::findSimilarGenePairs0(
    const string& geneSetName,
    const string& cellSetName,
    const string& similarGenePairsName,
    size_t k,                   // The maximum number of similar genes pairs to be stored for each gene.
    double similarityThreshold
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarGenePairs0 begins." << endl;

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }
    const GeneId geneCount = geneSet.size();

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Create the expression matrix subset for this gene set and cell set.
    cout << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset-" + similarGenePairsName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);



    // Create a dense expression vector for each gene.
    // All indices are local to the gene set and cell set.
    cout << timestamp << "Creating dense expression vectors." << endl;
    vector< vector<float> > v(geneCount, vector<float>(cellCount, 0.));
    for(CellId cellId=0; cellId!=cellCount; ++cellId) {
        for(const auto& p: expressionMatrixSubset.cellExpressionCounts[cellId]) {
            const GeneId geneId = p.first;
            const float count = p.second;
            v[geneId][cellId] = count;
        }
    }

    // Shift and normalize the dense expression vector of each gene
    // to zero mean and unit variance. This way regression coefficients can
    // be computed as simple scalar products.
    for(GeneId geneId=0; geneId!=geneCount; geneId++) {
        vector<float>& x = v[geneId];
        double sum = 0.;
        for(float count: x) {
            sum += count;
        }
        const float average = float(sum / cellCount);
        for(float& count: x) {
            count -= average;
        }
        double sum2 = 0.;
        for(float count: x) {
            sum2 += count * count;
        }
        const float factor = float(1./sqrt(sum2));
        for(float& count: x) {
            count *= factor;
        }
    }


    // Vectors to store the similar genes for each gene.
    vector< vector< pair<GeneId, float> > > similarGenes(geneCount);



    // Loop over gene pairs.
    size_t pairsDone = 0;
    size_t pairsStored = 0.;
    const size_t totalPairCount = (size_t(geneCount) * (size_t(geneCount) - 1)) / 2;
    cout << timestamp << "Loop over gene pairs begins." << endl;
    const size_t messageFrequency = size_t(1.e10/double(cellCount));
    for(GeneId geneId0=1; geneId0!=geneCount; geneId0++) {
        vector<float>& x0 = v[geneId0];
        for(GeneId geneId1=0; geneId1!=geneId0; geneId1++, ++pairsDone) {
            if(pairsDone>0 && (pairsDone % messageFrequency) == 0) {
                cout << timestamp << 100.*double(pairsDone)/double(totalPairCount) << "% done. ";
                cout << 100.*double(pairsStored)/double(pairsDone) << "% of gene pairs were stored.";
                cout << endl;
            }
            vector<float>& x1 = v[geneId1];
            const float r = std::inner_product(x0.begin(), x0.end(), x1.begin(), 0.f);
            if(r > similarityThreshold) {
                ++pairsStored;
                similarGenes[geneId0].push_back(make_pair(geneId1, r));
                similarGenes[geneId1].push_back(make_pair(geneId0, r));
            }
        }
    }



    // For each gene, keep only the k best and sort them.
    for(vector< pair<GeneId, float> >& v: similarGenes) {
        keepBest(v, k, OrderPairsBySecondGreater< pair<GeneId, float> >());
        sort(v.begin(), v.end(), OrderPairsBySecondGreater< pair<GeneId, float> >());
    }




    cout << timestamp << "ExpressionMatrix::findSimilarGenePairs0 ends." << endl;
}
