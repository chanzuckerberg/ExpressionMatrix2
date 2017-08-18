#include "ExpressionMatrix.hpp"
#include "ExpressionMatrixSubset.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "fstream.hpp"



// Compute a histogram of the difference between approximate and exact similarity,
// looping over all pairs. This is O(N**2) slow.
void ExpressionMatrix::analyzeAllPairs() const
{
    // Vector to contain the histogram of the difference between approximate and exact similarity.
    // The difference will always be in the range -2 to 2.
    const size_t binCount = 401;
    const double binWidth = 4. / double(binCount-1);
    const double binWidthInverse = 1. / binWidth;
    vector<uint64_t> histogram(binCount, 0ULL);

    // Loop over all pairs.
    for(CellId cellId0=0; cellId0!=cellCount()-1; cellId0++) {
        if((cellId0%100) == 0) {
            cout << timestamp << "Working on cell " << cellId0 << " of " << cells.size() << endl;
        }
        for(CellId cellId1=cellId0+1; cellId1!=cellCount(); cellId1++) {
            const double exactSimilarity = computeCellSimilarity(cellId0, cellId1);
            const double approximateSimilarity = computeApproximateCellSimilarity(cellId0, cellId1);
            const double delta = approximateSimilarity - exactSimilarity;
            size_t iDelta = std::lrint((delta+2.) * binWidthInverse);
            iDelta = max(iDelta, 0UL);
            iDelta = min(iDelta, binCount);
            ++(histogram[iDelta]);
        }
    }


    // Write out the histogram.
    ofstream csvOut("AllPairsHistogram.csv");
    for(size_t iDelta=0; iDelta<histogram.size(); iDelta++) {
        const double delta = double(iDelta) * binWidth - 2.;
        csvOut << iDelta << ",";
        csvOut << delta << ",";
        csvOut << histogram[iDelta] << "\n";
    }

}



// Find similar cell pairs by looping over all pairs.
// This can use an exact or approximate mode of operation, depending
// on the value of the last argument.
// In the exact mode of operation, all gene expression counts for both cells
// (stored in cellExpressionCounts) are taken into account.
// In the approximate mode of operation, only the largest gene expression counts for both cells
// (stored in largeCellExpressionCounts) are taken into account.
// The number of largest gene expression counts to be stored for each cell
// was specified when each cell was added, as an argument to addCell or addCells.
// Both the exact an the approximate mode of operation are O(N**2) slow
// because they loop over cell pairs. However the approximate mode of operation
// is usually much faster than the exact mode of operation
// (typically 10 times faster or better when the number of stored
// largest expression counts for each cell is 100.
void ExpressionMatrix::findSimilarPairs0(
    const string& cellSetName,	// The name of the cell set to be used.
    const string& name,         // The name of the SimilarPairs object to be created.
    size_t k,                   // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold, // The minimum similarity for a pair to be stored.
    bool useExactSimilarity     // Use exact of approximate cell similarity computation.
    )
{
    if(cellCount() == 0) {
        cout << "There are no cells. Skipping findSimilarPairs0." << endl;
        return;
    }

    // Create the SimilarPairs object where we will store the pairs.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + name, k, *(it->second));

    // Loop over all pairs.
    for(CellId localCellId0=0; localCellId0!=similarPairs.cellCount()-1; localCellId0++) {
        if(localCellId0>0 && ((localCellId0%100) == 0)) {
            cout << timestamp << "Working on cell " << localCellId0 << " of " << similarPairs.cellCount() << endl;
        }
        const CellId globalCellId0 = similarPairs.getGlobalCellId(localCellId0);

        // Find all cells with similarity better than the specified threshold.
        for(CellId localCellId1=localCellId0+1; localCellId1!=similarPairs.cellCount(); localCellId1++) {
            const CellId globalCellId1 = similarPairs.getGlobalCellId(localCellId1);
            const double similarity =
                useExactSimilarity ? computeCellSimilarity(globalCellId0, globalCellId1) : computeApproximateCellSimilarity(globalCellId0, globalCellId1);

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



// Version of findSimilarPairs0 that uses a gene set to find pairs of similar cells.
// Contrary to findSimilarPairs0, this only supports the exact mode.
void ExpressionMatrix::findSimilarPairs0WithGeneSet(
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
