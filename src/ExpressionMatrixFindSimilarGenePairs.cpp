#include "ExpressionMatrix.hpp"
#include "ExpressionMatrixSubset.hpp"
#include "heap.hpp"
#include "SimilarGenePairs.hpp"
#include "timestamp.hpp"
#include "tokenize.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "fstream.hpp"
#include <numeric>




void ExpressionMatrix::findSimilarGenePairs0(
    const string& geneSetName,
    const string& cellSetName,
    NormalizationMethod normalizationMethod,
    const string& similarGenePairsName,
    size_t k,                   // The maximum number of similar genes pairs to be stored for each gene.
    double similarityThreshold,
    bool writeCsv
    )
{
    findSimilarGenePairs0(cout, geneSetName, cellSetName,
        normalizationMethod, similarGenePairsName, k, similarityThreshold,
        writeCsv);
}
void ExpressionMatrix::findSimilarGenePairs0(
    ostream& s,
    const string& geneSetName,
    const string& cellSetName,
    NormalizationMethod normalizationMethod,
    const string& similarGenePairsName,
    size_t k,                   // The maximum number of similar genes pairs to be stored for each gene.
    double similarityThreshold,
    bool writeCsv
    )
{
    s << timestamp << "ExpressionMatrix::findSimilarGenePairs0 begins." << endl;
    s << "Gene set: " << geneSetName << endl;
    s << "Cell set: " << cellSetName << endl;
    s << "Normalization method: " << normalizationMethodToLongString(normalizationMethod) << endl;

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
    s << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset-" + similarGenePairsName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);



    // Create a dense expression vector for each gene.
    // All indices are local to the gene set and cell set.
    s << timestamp << "Creating dense expression vectors." << endl;
    vector< vector<float> > v;
    expressionMatrixSubset.getDenseRepresentation(v, normalizationMethod);



#if 0
    // The code under ifdef was moved to ExpressionMatrixSubset::getDenseRepresentation.
    vector< vector<float> > v(geneCount, vector<float>(cellCount, 0.));
    for(CellId cellId=0; cellId!=cellCount; ++cellId) {
        for(const auto& p: expressionMatrixSubset.cellExpressionCounts[cellId]) {
            const GeneId geneId = p.first;
            const float count = p.second;
            v[geneId][cellId] = count;
        }
    }



    // Normalize the expression vector of each cell, if requested.
    if(normalizationMethod != NormalizationMethod::none) {
        CZI_ASSERT(normalizationMethod != NormalizationMethod::Invalid);
        for(CellId cellId=0; cellId!=cellCount; cellId++) {
            const float factor =
                (normalizationMethod==NormalizationMethod::L1) ?
                    float(1. / expressionMatrixSubset.sums[cellId].sum1) :
                    float(1. / sqrt(expressionMatrixSubset.sums[cellId].sum2));
            for(GeneId geneId=0; geneId!=geneCount; geneId++) {
                v[geneId][cellId] *= factor;
            }
        }
    }
#endif


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

    // Open the csv file, if requested.
    ofstream csv;
    if(writeCsv) {
        csv.open(similarGenePairsName + ".csv");
    }

    // Loop over gene pairs.
    size_t pairsDone = 0;
    size_t pairsStored = 0.;
    const size_t totalPairCount = (size_t(geneCount) * (size_t(geneCount) - 1)) / 2;
    s << timestamp << "Loop over gene pairs begins." << endl;
    const size_t messageFrequency = size_t(1.e10/double(cellCount));
    for(GeneId geneId0=1; geneId0!=geneCount; geneId0++) {
        vector<float>& x0 = v[geneId0];
        for(GeneId geneId1=0; geneId1!=geneId0; geneId1++, ++pairsDone) {
            if(pairsDone>0 && (pairsDone % messageFrequency) == 0) {
                s << timestamp << 100.*double(pairsDone)/double(totalPairCount) << "% done. ";
                s << 100.*double(pairsStored)/double(pairsDone) << "% of gene pairs were stored.";
                s << endl;
            }
            vector<float>& x1 = v[geneId1];
            const float r = std::inner_product(x0.begin(), x0.end(), x1.begin(), 0.f);
            if(r > similarityThreshold) {
                ++pairsStored;
                similarGenes[geneId0].push_back(make_pair(geneId1, r));
                similarGenes[geneId1].push_back(make_pair(geneId0, r));
            }
            if(writeCsv) {
                csv << geneNames[geneId0] << ",";
                csv << geneNames[geneId1] << ",";
                csv << r << "\n";
                csv << geneNames[geneId1] << ",";
                csv << geneNames[geneId0] << ",";
                csv << r << "\n";
            }
        }
    }



    // For each gene, keep only the k best and sort them.
    s << timestamp << "Keeping best " << k << " pairs for each gene." << endl;
    size_t totalKept = 0;
    for(vector< pair<GeneId, float> >& v: similarGenes) {
        keepBest(v, k, OrderPairsBySecondGreater< pair<GeneId, float> >());
        sort(v.begin(), v.end(), OrderPairsBySecondGreater< pair<GeneId, float> >());
        totalKept += v.size();
    }
    s << "Average number of pairs kept per gene is " << double(totalKept)/geneCount << endl;


    // Create the SimilarGenePairs object.
    s << timestamp << "Permanently storing the similar gene pairs." << endl;
    SimilarGenePairs similarGenePairs(directoryName, similarGenePairsName,
        geneSetName, cellSetName, k, normalizationMethod, similarGenes);

    s << timestamp << "ExpressionMatrix::findSimilarGenePairs0 ends." << endl;
}



// Get a list of the currently available sets of similar gene pairs.
void ExpressionMatrix::getAvailableSimilarGenePairs(
    vector<string>& availableSimilarGenePairs) const
{

    const string fileNamePrefix = directoryName + "/SimilarGenePairs-";
    const string fileNameSuffix = "-Info";
    const vector<string> directoryContents = filesystem::directoryContents(directoryName);
    for(string name: directoryContents) {
        // Here, name contains the entire file name.
        if(stripPrefixAndSuffix(fileNamePrefix, fileNameSuffix, name)) {
            // Here, now contains just the similar pairs set name.
            availableSimilarGenePairs.push_back(name);
        }
    }

}


// Remove a similar gene pairs object given its name.
// This throws an exception if the requested SimilarGenePairs object does not exist.
void ExpressionMatrix::removeSimilarGenePairs(const string& name)
{
    try {
        SimilarGenePairs similarGenePairs(directoryName, name, false);
        similarGenePairs.remove();
    } catch(runtime_error e) {
        cout << e.what() << endl;
        throw runtime_error("Error removing similar gene pairs object " + name);
    }
}

