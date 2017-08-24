#ifndef CZI_EXPRESSION_MATRIX2_SIMILAR_PAIRS_HPP
#define CZI_EXPRESSION_MATRIX2_SIMILAR_PAIRS_HPP

#include "CellSets.hpp"
#include "GeneSet.hpp"
#include "Ids.hpp"
#include "MemoryAsContainer.hpp"
#include "MemoryMappedObject.hpp"
#include "MemoryMappedVector.hpp"

#include "string.hpp"
#include "utility.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class SimilarPairs;
    }
}



// A class to store similar pairs of cells.
// For each cell, we only store up to k similar cells.
// Therefore, if (cellId0, cellId1) is a stored pair,
// (cellId1, cellId0) is not guaranteed to also be a stored pair.

// IMPORTANT WARNING ABOUT CELL IDS:
// Note that all cell ids used and stored by this class are not global cell ids
// as seen by the ExpressionMatrix object. They are local cell ids to the cell set
// vector used by the SimilarPairs object.

// Only cell pairs where both cells are in the specified cell set are stored.

class ChanZuckerberg::ExpressionMatrix2::SimilarPairs {
public:

    // Create a new SimilarPairs object.
    SimilarPairs(
        const string& name,
        size_t k,
        const GeneSet& geneSet,
        const CellSet& cellSet);

    // Access an existing SimilarPairs object.
    SimilarPairs(const string& name);

    // The type used to store a cell similarity.
    typedef float CellSimilarity;

    // The type used to store a single pair. (The first cell is implied).
    typedef pair<CellId, CellSimilarity> Pair;

    // Returns a MemoryAsContainer containing the pairs for a given cell.
    // This can be used to easily iterate over the pairs for a given cell.
    MemoryAsContainer<const Pair>operator[](CellId cellId) const
    {
        const Pair* begin = similarPairs.begin() + cellId * k();
        const Pair* end = begin + usedCount[cellId];
        return  MemoryAsContainer<const Pair>(begin, end);
    }

    Pair* begin(CellId cellId)
    {
        return similarPairs.begin() + cellId * k();
    }
    const Pair* begin(CellId cellId) const
    {
        return similarPairs.begin() + cellId * k();
    }

    Pair* end(CellId cellId)
    {
        return begin(cellId) + size(cellId);
    }
    const Pair* end(CellId cellId) const
    {
        return begin(cellId) + size(cellId);
    }

    size_t size(CellId cellId) const
    {
        return usedCount[cellId];
    }

    // Return k, the maximum number of pairs stored for each cell.
    size_t k() const
    {
        return info->k;
    }

    // Add a pair.
    // This might or might not be stored, depending on the number
    // of pairs already stored for cellId0 and cellId1.
    // If could be stored or not for one or both of cellId0 and cellId1.
    void add(CellId cellId0, CellId cellId1, double similarity);

    // Return true if CellId1 is currently listed among the pairs
    // similar to CellId0.
    // Note that this function is not symmetric under a swap of cellId0 and cellid1.
    bool exists(CellId cellId0, CellId cellId1) const;

    // Sort the similar pairs for each cell by decreasing similarity.
    void sort();

    // Given a cell id local to this SimilarPairs object
    // (index in the cell set vector), return the corresponding
    // global CellId (index in the cells vector of the ExpressionMatrix object).
    // See the warning about cell ids before the definition of this class.
    CellId getGlobalCellId(CellId localCellId) const
    {
    	CZI_ASSERT(localCellId < cellSet.size());
    	return cellSet[localCellId];
    }

    // Inverse of the above.
    // Returns the local cell id corresponding to a global cell id.
    // Returns invalidCellId if none found.
    CellId getLocalCellId(CellId globalCellId) const
    {
    	const auto it = std::lower_bound(cellSet.begin(), cellSet.end(), globalCellId);
    	if(it == cellSet.end() || *it != globalCellId) {
    		return invalidCellId;
    	} else {
    		return CellId(it - cellSet.begin());
    	}
    }

    CellId cellCount() const
    {
        return CellId(cellSet.size());
    }

    const GeneSet& getGeneSet() const
    {
        return geneSet;
    }
    const CellSet& getCellSet() const
    {
        return cellSet;
    }

private:

    // All the pairs we could possibly store (k of them for each cell).
    // Stored contiguously, with the first k referring to cell0,
    // and so on.
    MemoryMapped::Vector<Pair> similarPairs;

    // Vector to hold the number of pairs actually stored for each cell.
    MemoryMapped::Vector<uint32_t> usedCount;

    // We also store, for each cell, a pair containing:
    // - The index, relative to begin() for this cell, of the pair with lowest similarity
    //   currently stored for this cell.
    // - The value of the lowest similarity currently stored for this cell.
    // If no pairs are stored for this cell, this stores a pair containing
    // (std::numeric_limits<uint32_t>::max() std::numeric_limits<CellSimilarity>::max()).
    MemoryMapped::Vector< pair<uint32_t, CellSimilarity> > lowestStoredSimilarityInfo;

    // Small size information about this table of similar pairs is stored
    // in a memory mapped object.

    class Info {
    public:
        // The maximum number of similar cells stored for each cell.
        size_t k;

        // The number of cells for which this table was constructed.
        CellId cellCount;
    };
    MemoryMapped::Object<Info> info;

    // The gene set and cell set used by this SimilarPairs object.
    // Tese are copies of the gene set and cell set passed to the constructor
    // when the SimilarPairs object was created.
    // These copies are owned by the SimilarPairs object.
    GeneSet geneSet;
    CellSet cellSet;

    // Add a pair (low level version).
    void add(CellId, Pair);
};


#endif
