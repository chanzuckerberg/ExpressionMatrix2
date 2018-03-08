#ifndef CZI_EXPRESSION_MATRIX2_SIMILAR_GENE_PAIRS_HPP
#define CZI_EXPRESSION_MATRIX2_SIMILAR_GENE_PAIRS_HPP

#include "CellSets.hpp"
#include "GeneSet.hpp"
#include "Ids.hpp"
#include "MemoryAsContainer.hpp"
#include "MemoryMappedObject.hpp"
#include "MemoryMappedVector.hpp"
#include "ShortStaticString.hpp"

#include "string.hpp"
#include "utility.hpp"


namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class SimilarGenePairs;
    }
}



// A class to store similar pairs of genes.
// For each cell, we only store up to k similar cells.
// Therefore, if (cellId0, cellId1) is a stored pair,
// (cellId1, cellId0) is not guaranteed to also be a stored pair.

// IMPORTANT WARNING ABOUT GENE IDS:
// Note that all gene ids used and stored by this class are not global gene ids
// as seen by the ExpressionMatrix object. They are local gene ids to the gene set
// vector used by the SimilarGenePairs object.

// Only gene pairs where both genes are in the specified gene set are stored.

class ChanZuckerberg::ExpressionMatrix2::SimilarGenePairs {
public:

    // The type used to store a cell similarity.
    using GeneSimilarity = float;

    // The type used to store a single pair. (The first gene is implied).
    using Pair = pair<GeneId, GeneSimilarity>;


    // Create a new SimilarGenePairs object.
    SimilarGenePairs(
        const string& directoryName,
        const string& similarGenePairsName,
        const string& geneSetName,
        const string& cellSetName,
        size_t k,
        const vector< vector<Pair> >&);

    // Access an existing SimilarGenePairs object.
    SimilarGenePairs(
        const string& directoryName,
        const string& similarGenePairsName,
        bool allowReadOnly);

    // Returns a MemoryAsContainer containing the pairs for a given gene.
    // This can be used to easily iterate over the pairs for a given gene.
    MemoryAsContainer<const Pair>operator[](GeneId geneId) const
    {
        const Pair* begin = pairs.begin() + geneId * k();
        const Pair* end = begin + geneInfo[geneId].usedCount;
        return  MemoryAsContainer<const Pair>(begin, end);
    }

    Pair* begin(GeneId geneId)
    {
        return pairs.begin() + geneId * k();
    }
    const Pair* begin(GeneId geneId) const
    {
        return pairs.begin() + geneId * k();
    }

    Pair* end(GeneId geneId)
    {
        return begin(geneId) + size(geneId);
    }
    const Pair* end(GeneId geneId) const
    {
        return begin(geneId) + size(geneId);
    }

    size_t size(GeneId geneId) const
    {
        return geneInfo[geneId].usedCount;
    }

    // Return k, the maximum number of pairs stored for each cell.
    size_t k() const
    {
        return info->k;
    }


    GeneId geneCount() const
    {
        return GeneId(geneSet.size());
    }

    const GeneSet& getGeneSet() const
    {
        return geneSet;
    }

    void remove();

private:

    // All the pairs we could possibly store (k of them for each gene).
    // Stored contiguously, with the first k referring to gene0,
    // and so on.
    MemoryMapped::Vector<Pair> pairs;



    // Some information stored for each gene.
    class GeneInfo {
    public:

        // The number of pairs stored so far for this gene (can be up to k).
        uint32_t usedCount;

    };
    MemoryMapped::Vector<GeneInfo> geneInfo;




    // Small size information about this table of similar pairs is stored
    // in a memory mapped object.
    class Info {
    public:
        // The maximum number of similar genes stored for each cell.
        size_t k;

        // The name and hash of the gene set used by this SimilarGenePairs object.
        StaticString255 geneSetName;
        uint64_t geneSetHash;

        // The name and hash of the cell set used by this SimilarGenePairs object.
        StaticString255 cellSetName;
        uint64_t cellSetHash;
    };
    MemoryMapped::Object<Info> info;

    static string getPathBaseName(
        const string& directoryName,
        const string& similarPairsName
        );


    // The gene set and cell set used by this SimilarGenePairs object.
    GeneSet geneSet;
    CellSet cellSet;
    void accessGeneSet(
        const string& directoryName,
        const string& geneSetName);
    void accessCellSet(
        const string& directoryName,
        const string& cellSetName);

};


#endif
