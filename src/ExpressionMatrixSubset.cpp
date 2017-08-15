// Class ExpressionMatrixSubset is used to store expression counts for a
// subset of cells and a subset of genes.

#include "ExpressionMatrixSubset.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;


ExpressionMatrixSubset::ExpressionMatrixSubset(
    const string& directoryName,
    const GeneSet& geneSet,
    const CellSet& cellSet,
    const CellExpressionCounts& globalExpressionCounts) :
        geneSet(geneSet), cellSet(cellSet)
{
    // Sanity checks.
    CZI_ASSERT(boost::filesystem::is_directory(directoryName));
    CZI_ASSERT(std::is_sorted(geneSet.begin(), geneSet.end()));
    CZI_ASSERT(std::is_sorted(cellSet.begin(), cellSet.end()));

    // Initialize the cell expression counts for the ExpressionMatrixSubset.
    cellExpressionCounts.createNew(directoryName);

    // Loop over cells in the subset.
#if 0
    for(CellId localCellId=0; localCellId!=cellSet.size(); localCellId++) {
        const CellId globalCellId = cellSet[localCellId];
        // Loop over all expression counts for this cell.
        for(const auto& p: globalExpressionCounts[globalCellId]) {
            const GeneId globalGeneId = p.first;
            const float count = p.second;
            CZI_ASSERT(0);
        }
    }
#endif

}

