#include "SimilarPairs.hpp"
#include "algorithm.hpp"
#include "orderPairs.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



// Create a new SimilarPairs object.
SimilarPairs::SimilarPairs(
    const string& name,
    size_t k,
    CellId cellCount)
{
    info.createNew(name + "-Info");
    info->k = k;
    info->cellCount = cellCount;

    similarPairs.createNew(name + "-Pairs", k*cellCount);
    usedCount.createNew(name + "-UsedCounts", k*cellCount);
    fill(usedCount.begin(), usedCount.end(), 0);

}



// Access an existing SimilarPairs object.
SimilarPairs::SimilarPairs(const string& name)
{
    info.accessExistingReadOnly(name + "-Info");
    similarPairs.accessExistingReadOnly(name + "-Pairs");
    usedCount.accessExistingReadOnly(name + "-UsedCounts");
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
    if(n < k()) {

        // There are unused slots. Just check if this pair already exists.
        for(uint32_t i=0; i<n; i++) {
            if(pairs[i].first == otherCellId) {
               return;  // Already exists.
            }
        }
        // Add it to one of the unused slots.
        pairs[n++] = pair;
        return;

    } else {

        // There are no unused slots. Check if this pair already exists,
        // and at the same time find the slot with the lowest similarity.
        uint32_t iLowestSimilarity = 0;
        double lowestSimilarity = std::numeric_limits<double>::max();
        for(uint32_t i=0; i<n; i++) {
            const Pair& existingPair = pairs[i];
            if(existingPair.first == otherCellId) {
               return;  // Already exists.
            }
            const double similarity = existingPair.second;
            if(similarity < lowestSimilarity) {
                iLowestSimilarity = i;
                lowestSimilarity = similarity;
            }
        }
        // Replace the entry in the slot with lowest similarity.
        if(pair.second > lowestSimilarity) {
            pairs[iLowestSimilarity] = pair;
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

