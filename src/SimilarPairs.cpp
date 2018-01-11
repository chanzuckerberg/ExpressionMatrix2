#include "SimilarPairs.hpp"
#include "algorithm.hpp"
#include "heap.hpp"
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

    similarPairs.createNew(name + "-Pairs", k*size_t(cellCount));

    // Initialize the cellInfo vector.
    cellInfo.createNew(name + "-CellInfo", cellCount);
    for(CellInfo& info: cellInfo) {
        info.usedCount = 0;
        info.lowestSimilarityIndex = std::numeric_limits<uint32_t>::max();
        info.lowestSimilarity = std::numeric_limits<CellSimilarity>::max();
    }

    // Make copies of the gene set and cell set. The copies are owned by the SimilarPairs object.
    geneSetArgument.makeCopy(geneSet, name + "-GeneSet");
    cellSetArgument.makeCopy(cellSet, name + "-CellSet");

}



// Access an existing SimilarPairs object.
SimilarPairs::SimilarPairs(const string& name, bool allowReadOnly)
{
    info.accessExistingReadOnly(name + "-Info");
    similarPairs.accessExistingReadOnly(name + "-Pairs");
    cellInfo.accessExistingReadOnly(name + "-CellInfo");
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
void SimilarPairs::addNoDuplicateCheck(CellId cellId0, CellId cellId1, double similarity)
{
    addNoDuplicateCheck(cellId0, make_pair(cellId1, similarity));
    addNoDuplicateCheck(cellId1, make_pair(cellId0, similarity));
}
void SimilarPairs::addNoDuplicateCheckUsingHeap(CellId cellId0, CellId cellId1, double similarity)
{
    addNoDuplicateCheckUsingHeap(cellId0, make_pair(cellId1, similarity));
    addNoDuplicateCheckUsingHeap(cellId1, make_pair(cellId0, similarity));
}



// This only adds to the list for cellId0.
void SimilarPairs::addUnsymmetric(CellId cellId0, CellId cellId1, double similarity)
{
    add(cellId0, make_pair(cellId1, similarity));
}
void SimilarPairs::addUnsymmetricNoDuplicateCheck(CellId cellId0, CellId cellId1, double similarity)
{
    addNoDuplicateCheck(cellId0, make_pair(cellId1, similarity));
}
void SimilarPairs::addUnsymmetricNoDuplicateCheckUsingHeap(CellId cellId0, CellId cellId1, double similarity)
{
    addNoDuplicateCheckUsingHeap(cellId0, make_pair(cellId1, similarity));
}


// Low level version of function to add a pair.
// The pair or might not be stored, depending on the number
// of pairs already stored for the given cell and whether it already exists.
// We could structure storage as a heap to reduce the insertion cost to O(logM),
// but there is no point because we need to do a linear search anyway
// to check for existence.
void SimilarPairs::add(CellId cellId, Pair pair)
{
    CellInfo& info = cellInfo[cellId];
    const CellId otherCellId = pair.first;
    const uint32_t n = info.usedCount;
    Pair* pairs = begin(cellId);
    if(n < k()) {

        // There are unused slots. Just check if this pair already exists.
        for(uint32_t i=0; i<n; i++) {
            if(pairs[i].first == otherCellId) {
               return;  // Already exists.
            }
        }
        // Update the lowest stored similarity info for this cell.
        if(pair.second < info.lowestSimilarity) {
            info.lowestSimilarityIndex = n;
            info.lowestSimilarity = pair.second;
        }

        // Add it to one of the unused slots.
        pairs[n] = pair;
        ++info.usedCount;

        return;

    } else {

        // There are no unused slots.

        // If the similarity of this pair is less than the lowest stored
        // similarity for this cell, do nothing.
        // This way we avoid a scan of the stored pairs for this cell.
        if(pair.second <= info.lowestSimilarity) {
            return;
        }


        // Check if this pair already exists.
        for(uint32_t i = 0; i < n; i++) {
            const Pair& existingPair = pairs[i];
            if(existingPair.first == otherCellId) {
                return;  // Already exists.
            }
        }

        // Store the pair in the slot containing the pair with the lowest similarity.
        pairs[info.lowestSimilarityIndex] = pair;

        // Update the lowest similarity.
        info.lowestSimilarityIndex = std::numeric_limits<uint32_t>::max();
        info.lowestSimilarity = std::numeric_limits<CellSimilarity>::max();
        for(uint32_t i=0; i<n; i++) {
            const Pair& existingPair = pairs[i];
            if(existingPair.second < info.lowestSimilarity) {
                info.lowestSimilarityIndex = i;
                info.lowestSimilarity = existingPair.second;
            }
        }

        return;
    }
}



// This version does a linear search to update the position of lowest
// similarity item.
void SimilarPairs::addNoDuplicateCheck(CellId cellId, Pair pair)
{
    CellInfo& info = cellInfo[cellId];
    const uint32_t n = info.usedCount;
    if(n < k()) {

        // Update the lowest stored similarity info for this cell.
        if(pair.second < info.lowestSimilarity) {
            info.lowestSimilarityIndex = n;
            info.lowestSimilarity = pair.second;
        }

        // Add it to one of the unused slots.
        Pair* pairs = begin(cellId);
        pairs[n] = pair;
        ++info.usedCount;

        return;

    } else {

        // There are no unused slots.

        // If the similarity of this pair is less than the lowest stored
        // similarity for this cell, do nothing.
        // This way we avoid a scan of the stored pairs for this cell.
        if(pair.second <= info.lowestSimilarity) {
            return;
        }


        // Store the pair in the slot containing the pair with the lowest similarity.
        Pair* pairs = begin(cellId);
        pairs[info.lowestSimilarityIndex] = pair;

        // Update the lowest similarity.
        info.lowestSimilarityIndex = std::numeric_limits<uint32_t>::max();
        info.lowestSimilarity = std::numeric_limits<CellSimilarity>::max();
        for(uint32_t i=0; i<n; i++) {
            const Pair& existingPair = pairs[i];
            if(existingPair.second < info.lowestSimilarity) {
                info.lowestSimilarityIndex = i;
                info.lowestSimilarity = existingPair.second;
            }
        }

        return;
    }
}



void SimilarPairs::addUnsymmetricNoCheck(CellId cellId0, CellId cellId1, double similarity)
{
    CellInfo& info0 = cellInfo[cellId0];
    uint32_t& n0 = info0.usedCount;
    CZI_ASSERT(n0 < k());

    Pair* pairs = begin(cellId0);
    Pair& pair = pairs[n0];
    pair.first = cellId1;
    pair.second = float(similarity);
    ++n0;

}


// Version that uses a heap to avoid linear searches.
void SimilarPairs::addNoDuplicateCheckUsingHeap(CellId cellId, Pair pair)
{
    CellInfo& info = cellInfo[cellId];
    const uint32_t n = info.usedCount;
    const size_t kk = k();
    if(n < kk) {

        // Store it in the next unused slot.
        Pair* pairs = begin(cellId);
        pairs[n] = pair;
        ++info.usedCount;

        // If we filled up all the slots, turn it into a heap,
        // so the next insertions will work.
        if(info.usedCount == kk) {
            std::make_heap(pairs, pairs+kk, OrderPairsBySecondGreater<Pair>());
        }

        return;

    } else {

#if 0
        cout << "Adding " << pair.first << " " << pair.second << endl;
        cout << "Heap is: " << endl;
        for(auto it=begin(cellId); it!=end(cellId); ++it) {
            cout << it->first << " " << it->second << endl;
        }
#endif

        // The pair with the lowest similarity is the first one in the heap.
        // If the pair passed as an argument has lower similarity than that,
        // don't store it.
        Pair* b = begin(cellId);
        if(pair.second <= b->second) {
            // cout << "Skipped." << endl;
            return;
        }

        // The pair passed as an argument has higher similarity
        // than the lowest similarity currently stored for this cell.
        // Update the heap.
        Pair* e = b + kk;
#if 1
        std::pop_heap(b, e, OrderPairsBySecondGreater<Pair>());
        *(e-1) = pair;  // Store in the last slot. push_head will move it as required.
        std::push_heap(b, e, OrderPairsBySecondGreater<Pair>());
#else
        popAndPushHeap(b, e, pair, OrderPairsBySecondGreater<Pair>());
#endif

#if 0
        cout << "After adding, heap is: " << endl;
        for(auto it=begin(cellId); it!=end(cellId); ++it) {
            cout << it->first << " " << it->second << endl;
        }
#endif
        return;
    }
}


// Copy the pairs from the argument.
void SimilarPairs::copy(const vector< vector<Pair> >& v)
{
    CZI_ASSERT(v.size() == size_t(cellCount()));
    for(CellId cellId=0; cellId<cellCount(); cellId++) {
        const vector<Pair>& x = v[cellId];
        CZI_ASSERT(x.size() <= k());
        std::copy(x.begin(), x.end(), begin(cellId));
        cellInfo[cellId].usedCount = CellId(x.size());
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
        std::sort(begin(cellId), end(cellId), OrderPairsBySecondGreaterThenByFirstLess<Pair>());
    }

}


