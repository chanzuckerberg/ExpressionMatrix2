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

	void createNew(const string& name);
	void accessExisting(const string& name);

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
	size_t size() const
	{
		return geneVector.size();
	}


	const GeneId* begin() const
	{
		return geneVector.begin();
	}

	const GeneId* end() const
	{
		return geneVector.end();
	}

	GeneId operator[](size_t i) const
	{
		return geneVector[i];
	}

	void remove()
	{
		geneVector.remove();
		isGeneInSet.remove();
	}

	void getSortedGenes(vector<GeneId>&);

private:

	// The ids of the genes in this set.
	MemoryMapped::Vector<GeneId> geneVector;

	// Flag that tells us whether the genes vector is currently sorted by GeneId.
	bool isSorted = false;

	// Sort the genes vector and set the isSorted flag.
	void sort();

	// Vector that stores true fore each gene in the gene set
	// and false otherwise. Indexed by the GeneId.
	MemoryMapped::Vector<bool> isGeneInSet;
};


#endif
