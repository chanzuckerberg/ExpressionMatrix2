#include "GeneSet.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "algorithm.hpp"



void GeneSet::createNew(const string& name)
{
	geneVector.createNew(name + "-Genes", 0);
	isGeneInSet.createNew(name + "-Flags", 0);
}



void GeneSet::accessExisting(const string& name)
{
	geneVector.accessExistingReadWrite(name + "-Genes");
	isGeneInSet.accessExistingReadWrite(name + "-Flags");
	sort();
}



// Add a gene to the set.
// We lazily add it to the end of the gene vector,
// postponing the sorting of the vector.
void GeneSet::addGene(GeneId geneId)
{
	// Update the gene vector and its sorted flag.
	geneVector.push_back(geneId);
	isSorted = false;

	// Update the isGeneInSet vector of boolean flags.
	if(geneId >= isGeneInSet.size()) {
		const size_t oldSize = isGeneInSet.size();
		isGeneInSet.resize(geneId+1);
		const size_t newSize = isGeneInSet.size();
		fill(isGeneInSet.begin()+oldSize, isGeneInSet.begin()+newSize, false);
	}
	isGeneInSet[geneId] = true;
}



// Return true if a given GeneId belongs to the set, false otherwise.
bool GeneSet::contains(GeneId geneId)
{
	if(geneId < isGeneInSet.size()) {
		return isGeneInSet[geneId];
	} else {
		return false;
	}
}



// Get a reference to the vector of gene ids for the genes in this set.
// When this function returns, the vector is guaranteed to be sorted by id.
const MemoryMapped::Vector<GeneId>& GeneSet::genes()
{
	sort();
	return geneVector;

}



// Sort the genes vector and set the isSorted flag, if necessary.
void GeneSet::sort()
{
	if(!isSorted) {
		std::sort(geneVector.begin(), geneVector.end());
		isSorted = true;
	}
}



// Get a vector of genes in the set, sorted by GeneId.
void GeneSet::getSortedGenes(vector<GeneId>& sortedGenes)
{
    sort();
    sortedGenes.clear();
    sortedGenes.resize(geneVector.size());
    copy(geneVector.begin(), geneVector.end(), sortedGenes.begin());
}
