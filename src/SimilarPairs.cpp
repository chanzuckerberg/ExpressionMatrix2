#include "SimilarPairs.hpp"
#include "algorithm.hpp"
#include "orderPairs.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



// Create a new SimilarPairs object.
SimilarPairs::SimilarPairs(
    const string& name,
    size_t k,
    const GeneSet& geneSetArgument,
    const CellSet& cellSetArgument)
{
    const CellId cellCount = CellId(cellSetArgument.size());

    info.createNew(name + "-Info");
    info->k = k;
    info->cellCount = cellCount;

    similarPairs.createNew(name + "-Pairs", k*cellSetArgument.size());

    // Initialize the usedCount for each cell to zero.
    usedCount.createNew(name + "-UsedCounts", cellCount);
    fill(usedCount.begin(), usedCount.end(), 0);

    // Initialized the lowestStoredSimilarity for each cell to its maximum possible value..
    lowestStoredSimilarityInfo.createNew(name + "-LowestStoredSimilarityInfo", cellCount);
    fill(lowestStoredSimilarityInfo.begin(), lowestStoredSimilarityInfo.end(),
        make_pair(std::numeric_limits<uint32_t>::max(), std::numeric_limits<CellSimilarity>::max()));

    // Make copies of the gene set and cell set. The copies are owned by the SimilarPairs object.
    geneSetArgument.makeCopy(geneSet, name + "-GeneSet");
    cellSetArgument.makeCopy(cellSet, name + "-CellSet");

}



// Access an existing SimilarPairs object.
SimilarPairs::SimilarPairs(const string& name, bool allowReadOnly)
{
    info.accessExistingReadOnly(name + "-Info");
    similarPairs.accessExistingReadOnly(name + "-Pairs");
    usedCount.accessExistingReadOnly(name + "-UsedCounts");
    lowestStoredSimilarityInfo.accessExistingReadOnly(name + "-LowestStoredSimilarityInfo");
    geneSet.accessExisting(name + "-GeneSet", allowReadOnly);
    cellSet.accessExistingReadOnly(name + "-CellSet");
}



// Add a pair.
// This might or might not be stored, depending on the number
// of pairs already stored for cellId0 and cellId1.
// If could be stored or not for one or both of cellId0 and cellId1.
void SimilarPairs::add(CellId cellId0, CellId cellId1, double similarity)
{
    add(cellId0, make_pair(cellId1, similarity));
    add(cellId1, make_pair(cellId0, similarity));
}



// Low level version of function to add a pair.
// The pair or might not be stored, depending on the number
// of pairs already stored for the given cell and whether it already exists.
// We could structure storage as a heap to reduce the insertion cost to O(logM),
// but there is no point because we need to do a linear search anyway
// to check for existence.
void SimilarPairs::add(CellId cellId, Pair pair)
{
    const CellId otherCellId = pair.first;
    uint32_t& n = usedCount[cellId];
    Pair* pairs = begin(cellId);
    auto& lowestInfo = lowestStoredSimilarityInfo[cellId];
    if(n < k()) {

        // There are unused slots. Just check if this pair already exists.
        for(uint32_t i=0; i<n; i++) {
            if(pairs[i].first == otherCellId) {
               return;  // Already exists.
            }
        }
        // Update the lowest stored similarity info for this cell.
        if(pair.second < lowestInfo.second) {
        	lowestInfo.first = n;
        	lowestInfo.second = pair.second;
        }

        // Add it to one of the unused slots.
        pairs[n++] = pair;

        return;

    } else {

        // There are no unused slots.

    	// If the similarity of this pair is less than the lowest stored
    	// similarity for this cell, do nothing.
    	// This way we avoid a scan of the stored pairs for this cell.
    	if(pair.second <= lowestInfo.second) {
    		return;
    	}


    	// Check if this pair already exists.
        for(uint32_t i=0; i<n; i++) {
            const Pair& existingPair = pairs[i];
            if(existingPair.first == otherCellId) {
               return;  // Already exists.
            }
        }

        // Store the pair in the slot containing the pair with the lowest similarity.
        pairs[lowestInfo.first] = pair;

        // Update the lowest similarity.
        lowestInfo.first = std::numeric_limits<uint32_t>::max();
        lowestInfo.second = std::numeric_limits<CellSimilarity>::max();
        for(uint32_t i=0; i<n; i++) {
            const Pair& existingPair = pairs[i];
            if(existingPair.second < lowestInfo.second) {
               lowestInfo.first = i;
               lowestInfo.second = existingPair.second;
            }
        }

        return;
    }
}



// Return true if CellId1 is currently listed among the pairs
// similar to CellId0.
// Note that this function is not symmetric under a swap of cellId0 and cellid1.
bool SimilarPairs::exists(CellId cellId0, CellId cellId1) const
{
    for(const Pair& pair: (*this)[cellId0]) {
        if(pair.first == cellId1) {
            return true;
        }
    }
    return false;
}



// Sort the similar pairs for each cell by decreasing similarity.
void SimilarPairs::sort()
{
    for(CellId cellId=0; cellId<info->cellCount; cellId++) {
        std::sort(begin(cellId), end(cellId), OrderPairsBySecondGreater<Pair>());
    }

}

