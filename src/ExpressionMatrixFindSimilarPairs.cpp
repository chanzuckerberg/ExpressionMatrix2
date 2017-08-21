#include "ExpressionMatrix.hpp"
#include "ExpressionMatrixSubset.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "fstream.hpp"



// Find similar cell pairs by looping over all pairs,
// taking into account only genes in the specified gene set.
// This is O(N**2) slow because it loops over cell pairs.
void ExpressionMatrix::findSimilarPairs0(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold
    )
{
    // Locate the gene set.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;

    // Locate the cell set.
    const auto itCellSet = cellSets.cellSets.find(cellSetName);
    if(itCellSet == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const CellSet& cellSet = *(itCellSet->second);

    if(cellCount() == 0) {
        cout << "There are no cells. Skipping findSimilarPairs0." << endl;
        return;
    }

    // Create the SimilarPairs object where we will store the pairs.
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, k, cellSet);

    // Create the expression matrix subset for this gene set and cell set.
    const string expressionMatrixSubsetName = directoryName + "/tmp-ExpressionMatrixSubset-" + similarPairsName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);


    // Loop over all pairs.
    for(CellId localCellId0=0; localCellId0!=similarPairs.cellCount()-1; localCellId0++) {
        if(localCellId0>0 && ((localCellId0%100) == 0)) {
            cout << timestamp << "Working on cell " << localCellId0 << " of " << cellSet.size() << endl;
        }

        // Find all cells with similarity better than the specified threshold.
        for(CellId localCellId1=localCellId0+1; localCellId1!=similarPairs.cellCount(); localCellId1++) {
            const double similarity = expressionMatrixSubset.computeCellSimilarity(localCellId0, localCellId1);

            // If the similarity is sufficient, pass it to the SimilarPairs container,
            // which will make the decision whether to store it, depending on the
            // number of pairs already stored for cellId0 and cellId1.
            if(similarity > similarityThreshold) {
                similarPairs.add(localCellId0, localCellId1, similarity);
            }
        }
    }


    // Sort the similar pairs for each cell by decreasing similarity.
    similarPairs.sort();
}



// Dump cell to csv file a set of similar cell pairs.
void ExpressionMatrix::writeSimilarPairs(const string& name) const
{
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + name);
    ofstream csvOut("SimilarPairs-" + name + ".csv");
    csvOut << "Cell0,Cell1,Computed,Exact\n";

    for(CellId cellId0=0; cellId0<cellCount(); cellId0++) {
        for(const auto& pairs0: similarPairs[cellId0]) {
            csvOut << cellId0 << ",";
            csvOut << pairs0.first << ",";
            csvOut << pairs0.second << ",";
            csvOut << computeCellSimilarity(cellId0, pairs0.first) << "\n";
        }
    }

}



// Unit test for class ExpressionMatrixSubset.
void ExpressionMatrix::testExpressionMatrixSubset(CellId cellId0, CellId cellId1) const
{
    cout << "Test ExpressionMatrixSubset on cell ids " << cellId0 << " " << cellId1 << endl;

    // Locate the gene set.
    const auto itGeneSet = geneSets.find("AllGenes");
    CZI_ASSERT(itGeneSet != geneSets.end());
    const GeneSet& geneSet = itGeneSet->second;

    // Locate the cell set.
    const auto itCellSet = cellSets.cellSets.find("AllCells");
    CZI_ASSERT(itCellSet != cellSets.cellSets.end());
    const CellSet& cellSet = *(itCellSet->second);

    const string expressionMatrixSubsetName = directoryName + "/tmp-ExpressionMatrixSubset-Test";
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);

    cout << "Exact " << computeCellSimilarity(cellId0, cellId1) << endl << endl;
    cout << "Subset " << expressionMatrixSubset.computeCellSimilarity(cellId0, cellId1) << endl;
}
