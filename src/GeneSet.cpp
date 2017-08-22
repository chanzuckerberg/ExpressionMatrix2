#include "GeneSet.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "algorithm.hpp"



void GeneSet::createNew(const string& name)
{
    globalGeneIdVector.createNew(name + "-GlobalIds", 0);
    localGeneIdVector.createNew(name + "-LocalIds", 0);
}



void GeneSet::accessExisting(const string& name)
{
    globalGeneIdVector.accessExistingReadWrite(name + "-GlobalIds");
    localGeneIdVector.accessExistingReadWrite(name + "-LocalIds");
    sort();
}



// Make a copy of this gene set.
void GeneSet::makeCopy(GeneSet& copy, const string& newName) const
{
    globalGeneIdVector.makeCopy(copy.globalGeneIdVector, newName + "-GlobalIds");
    localGeneIdVector.makeCopy(copy.localGeneIdVector, newName + "-LocalIds");
    copy.sort();
}



// Add a gene to the set.
// We lazily add it to the end of the globalGeneId gene vector,
// postponing the sorting of the vector.
void GeneSet::addGene(GeneId geneId)
{
    if(geneId >= localGeneIdVector.size()) {
        const auto oldSize = localGeneIdVector.size();
        localGeneIdVector.resize(geneId+1);
        fill(localGeneIdVector.begin()+oldSize, localGeneIdVector.end(), invalidGeneId);
    }
    localGeneIdVector[geneId] = GeneId(globalGeneIdVector.size());
    globalGeneIdVector.push_back(geneId);
    isSorted = false;
}



// Return true if a given global GeneId belongs to the set, false otherwise.
bool GeneSet::contains(GeneId globalGeneId)
{
    return getLocalGeneId(globalGeneId) != invalidGeneId;
}



// Get a reference to the vector of gene ids for the genes in this set.
// When this function returns, the vector is guaranteed to be sorted by id.
const MemoryMapped::Vector<GeneId>& GeneSet::genes()
{
    sort();
    return globalGeneIdVector;
}



// Sort the genes vector and set the isSorted flag, if necessary.
void GeneSet::sort()
{
    if(!isSorted) {
        std::sort(globalGeneIdVector.begin(), globalGeneIdVector.end());
        isSorted = true;

        // Update the localGeneId vector to reflect the new order.
        fill(localGeneIdVector.begin(), localGeneIdVector.end(), invalidGeneId);
        for(GeneId localGeneId=0; localGeneId!=globalGeneIdVector.size(); localGeneId++) {
            const GeneId globalGeneId = globalGeneIdVector[localGeneId];
            localGeneIdVector[globalGeneId] = localGeneId;
        }
    }
}



// Get a vector of genes in the set, sorted by GeneId.
void GeneSet::getSortedGenes(vector<GeneId>& sortedGenes)
{
    sort();
    sortedGenes.clear();
    sortedGenes.resize(globalGeneIdVector.size());
    copy(globalGeneIdVector.begin(), globalGeneIdVector.end(), sortedGenes.begin());
}


// The comparison operator requires the two gene sets being compared
// to be sorted. It will assert if this is not the case.
bool GeneSet::operator==(const GeneSet& that) const
{
    CZI_ASSERT(isSorted);
    CZI_ASSERT(that.isSorted);
    return globalGeneIdVector == that.globalGeneIdVector;
}

