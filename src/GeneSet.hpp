#ifndef CZI_EXPRESSION_MATRIX2_GENE_SET_HPP
#define CZI_EXPRESSION_MATRIX2_GENE_SET_HPP

#include "Ids.hpp"
#include "MemoryMappedVector.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class GeneSet;
    }
}



// Class to describe a subset of all genes stored in the ExpressionMatrix object.
// The total number of genes will stay small, so we can store the gene set
// in multiple ways to facilitate processing.
class ChanZuckerberg::ExpressionMatrix2::GeneSet {
public:

    // Create a new GeneSet.
    void createNew(const string& name);

    // Access a previously created GeneSet.
    void accessExisting(const string& name);

    // Make a copy of this gene set.
    void makeCopy(GeneSet& copy, const string& newName) const;

    // Add a gene to the set.
    // The caller is responsible to make sure that a
    // gene is not added more than once.
    void addGene(GeneId);

    // Return true if a given GeneId belongs to the set, false otherwise.
    bool contains(GeneId);

    // Get a reference to the vector of gene ids for the genes in this set.
    // When this function returns, the vector is guaranteed to be sorted by id.
    const MemoryMapped::Vector<GeneId>& genes();

    // Return the number of genes in the set.
    GeneId size() const
    {
        return GeneId(globalGeneIdVector.size());
    }


    const GeneId* begin() const
    {
        return globalGeneIdVector.begin();
    }

    const GeneId* end() const
    {
        return globalGeneIdVector.end();
    }

    // Return the global gene id corresponding to a given local GeneId.
    GeneId getGlobalGeneId(GeneId localGeneId) const
    {
        CZI_ASSERT(localGeneId < size());
        return globalGeneIdVector[localGeneId];
    }

    // Return the local gene id corresponding to a given global GeneId,
    // or invalidGeneId of the specified global GeneId is not in this GeneSet.
    GeneId getLocalGeneId(GeneId globalGeneId) const
    {
        if(globalGeneId < localGeneIdVector.size()) {
            return localGeneIdVector[globalGeneId];
        } else {
            return invalidGeneId;
        }
    }

    void remove()
    {
        globalGeneIdVector.remove();
        localGeneIdVector.remove();
    }

    void getSortedGenes(vector<GeneId>&);

    // The comparison operator requires the two gene sets being compared
    // to be sorted. It will assert if this is not the case.
    bool operator==(const GeneSet&) const;

    // Sort the genes vector and set the isSorted flag.
    void sort();

    // Assert that the gene set is sorted.
    void assertIsSorted() const
    {
        CZI_ASSERT(isSorted);
    }

    // Set the isSorted flag without checking that the ste is sorted.
    // The caller is responsible for guaranteeing that the set is sorted.
    void forceSorted()
    {
        isSorted = true;
    }

private:

    // The global GeneId's of the genes in this set.
    // Indexed by the local GeneId.
    MemoryMapped::Vector<GeneId> globalGeneIdVector;

    // Flag that tells us whether the genes vector is currently sorted by GeneId.
    bool isSorted = false;

    // Vector indexed by the global GeneId that contains
    // the local GeneId for each gene in the gene set,
    // or invalidGeneId if that global gene id is not in the set.
    // Note that this needs to be recomputed when the GeneSet
    // gets sorted.
    MemoryMapped::Vector<GeneId> localGeneIdVector;


};


#endif
