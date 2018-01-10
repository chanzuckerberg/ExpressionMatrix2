// Classes used to implement the Charikar algorithm
// in ExpressionMatrix::findSimilarPairs6.
// See M. Charikar, "Similarity Estimation Techniques from Rounding Algorithms", 2002,
// section "5. Approximate Nearest neighbor Search in Hamming Space.",
// and ExpressionMatrix::findSimilarPairs6 in ExpressionMatrixLsh.cpp
// for more information.

#ifndef CZI_EXPRESSION_MATRIX2_CHARIKAR_HPP
#define CZI_EXPRESSION_MATRIX2_CHARIKAR_HPP



#include "BitSet.hpp"
#include "Ids.hpp"

#include "algorithm.hpp"
#include <limits>
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        namespace Charikar {
            class PermutationData;
            class Pointer;
        }
    }
}


// The information stored for each of the bit permutations.
class ChanZuckerberg::ExpressionMatrix2::Charikar::PermutationData {
public:

    // The signatures, with the bits shuffled according
    // to the bit permutation this data structure refers to,
    // and sorted in lexicographical order.
    // We only store the most significant permutedBitCount bits
    // of each signature.
    BitSets signatures;

    // The cell ids, in corresponding order.
    // All cell ids are local to the cell set in use.
    vector<CellId> cellIds;

    // The constructor allocates the memory but does not fill it in.
    PermutationData(CellId cellCount, size_t permutedWordCount) :
        signatures(cellCount, permutedWordCount),
        cellIds(cellCount),
        cellPositions(cellCount)
    {
    }

    // Store the position in the signatures and cellIds vectors
    // corresponding to each CellId.
    vector<size_t> cellPositions;
    void computeCellPositions()
    {
        const CellId cellCount = CellId(cellIds.size());
        CZI_ASSERT(cellPositions.size() == cellCount);
        fill(cellPositions.begin(), cellPositions.end(),
            std::numeric_limits<CellId>::max());
        for(size_t i=0; i!=cellCount; i++) {
            const CellId cellId = cellIds[i];
            cellPositions[cellId] = i;
        }
    }
};



// The class used to store a Charikar pointer.
// To run the Charikar algorithm for a cell, we store a priority queue
// of objects of this type.
class ChanZuckerberg::ExpressionMatrix2::Charikar::Pointer {
public:

    // The permutation that this pointer refers to.
    size_t permutationId;

    // The index in the signatures and cellIds vectors.
    size_t index;

    // Whether this pointer moves forward or backward.
    bool movesForward;

    // The number of identical prefix bits between
    // the pointed to signature and the permuted signature
    // of the cell we are looking for.
    size_t prefixLength;

    // Order by prefix length.
    // The priority queue will have as its top element
    // the one with the greatest prefix length.
    bool operator<(const Pointer& that) const
    {
        return prefixLength<that.prefixLength;
    }

};
#endif
