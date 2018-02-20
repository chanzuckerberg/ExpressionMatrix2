// Class ExpressionMatrixSubset is used to store expression counts for a
// subset of cells and a subset of genes.

#include "ExpressionMatrixSubset.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;


ExpressionMatrixSubset::ExpressionMatrixSubset(
    const string& name,
    const GeneSet& geneSet,
    const CellSet& cellSet,
    const CellExpressionCounts& globalExpressionCounts) :
        geneSet(geneSet), cellSet(cellSet)
{
    // Sanity checks.
    CZI_ASSERT(std::is_sorted(geneSet.begin(), geneSet.end()));
    CZI_ASSERT(std::is_sorted(cellSet.begin(), cellSet.end()));

    // Initialize the cell expression counts for the ExpressionMatrixSubset.
    cellExpressionCounts.createNew(name);

    // Loop over cells in the subset.
    for(CellId localCellId=0; localCellId!=cellSet.size(); localCellId++) {
        const CellId globalCellId = cellSet[localCellId];
        cellExpressionCounts.appendVector();

        // Loop over all expression counts for this cell.
        for(const auto& p: globalExpressionCounts[globalCellId]) {
            const GeneId globalGeneId = p.first;
            const GeneId localGeneId = geneSet.getLocalGeneId(globalGeneId);
            if(localGeneId == invalidGeneId) {
                continue;   // This gene is not in the gene set
            }
            const float count = p.second;
            cellExpressionCounts.append(make_pair(localGeneId, count));
        }
    }

    // Compute sums and sums of squares of the expression counts for all cells.
    computeSums();
}



// Compute sums and sums of squares of the expression counts for all cells.
void ExpressionMatrixSubset::computeSums()
{
    sums.resize(cellSet.size());
    for(CellId localCellId=0; localCellId<cellSet.size(); localCellId++) {
        Sum& sum = sums[localCellId];
        for(const auto& p: cellExpressionCounts[localCellId]) {
            const float& count = p.second;
            sum.sum1 += count;
            sum.sum2 += count*count;
        }
    }
}



ExpressionMatrixSubset::~ExpressionMatrixSubset()
{
    remove();
}



// Close and remove the supporting files.
void ExpressionMatrixSubset::remove()
{
    cellExpressionCounts.remove();
}



// Compute the similarity between two cells, identified by their ids
// local to our cell set, and using the stored expression counts
// (which reflect only genes in our gene set).
// This is similar to ExpressionMatrix:;computeCellSimilarity,
// but it used the expression counts stored in the ExpressionMatrixSubset
// instead of the global expression counts stored by class ExpressionMatrix.
double ExpressionMatrixSubset::computeCellSimilarity(CellId localCellId0, CellId localCellId1) const
{
    // Compute the scalar product of the expression counts for the two cells.
    typedef pair<GeneId, float>const* Iterator;
    const Iterator begin0 = cellExpressionCounts.begin(localCellId0);
    const Iterator end0 = cellExpressionCounts.end(localCellId0);
    const Iterator begin1 = cellExpressionCounts.begin(localCellId1);
    const Iterator end1 = cellExpressionCounts.end(localCellId1);
    Iterator it0 = begin0;
    Iterator it1 = begin1;
    double scalarProduct = 0.;
    while((it0 != end0) && (it1 != end1)) {
        const GeneId localGeneId0 = it0->first;
        const GeneId localGeneId1 = it1->first;

        if(localGeneId0 < localGeneId1) {
            ++it0;
        } else if(localGeneId1 < localGeneId0) {
            ++it1;
        } else {
            scalarProduct += it0->second * it1->second;
            ++it0;
            ++it1;
        }
    }

    // Compute the correlation coefficient.
    // See, for example, https://en.wikipedia.org/wiki/Correlation_and_dependence
    const double n = double(geneSet.size());
    const Sum& cell0Sum = sums[localCellId0];
    const Sum& cell1Sum = sums[localCellId1];
    const double numerator = n*scalarProduct - cell0Sum.sum1*cell1Sum.sum1;
    const double denominator = sqrt(
            (n*cell0Sum.sum2 - cell0Sum.sum1*cell0Sum.sum1) *
            (n*cell1Sum.sum2 - cell1Sum.sum1*cell1Sum.sum1)
            );

#if 0
    cout << "ExpressionMatrixSubset::computeCellSimilarity" << endl;
    cout << "n " << n << endl;
    cout << "scalarProduct " << scalarProduct << endl;
    cout << "cell0Sum.sum1 " << cell0Sum.sum1 << endl;
    cout << "cell1Sum.sum1 " << cell1Sum.sum1 << endl;
    cout << "cell0Sum.sum2 " << cell0Sum.sum2 << endl;
    cout << "cell1Sum.sum2 " << cell1Sum.sum2 << endl;
    cout << "numerator " << numerator << endl;
    cout << "denominator " << denominator << endl;
#endif

    return numerator / denominator;
}



// Get a dense representation of the expression matrix subset.
// Indexed by [localGeneId][localCellId].
// Note this means that the counts for all cells and a given
// gene are contiguous.
// Cell expression vectors are normalized as requested.
void ExpressionMatrixSubset::getDenseRepresentation(
    vector< vector<float> >& v,
    NormalizationMethod normalizationMethod) const
{
    v.clear();
    v.resize(geneCount(), vector<float>(cellCount(), 0.));

    for(CellId cellId=0; cellId!=cellCount(); ++cellId) {
        for(const auto& p: cellExpressionCounts[cellId]) {
            const GeneId geneId = p.first;
            const float count = p.second;
            v[geneId][cellId] = count;
        }
    }

    // Normalize the expression vector of each cell, if requested.
    if(normalizationMethod != NormalizationMethod::none) {
        CZI_ASSERT(normalizationMethod != NormalizationMethod::Invalid);
        for(CellId cellId=0; cellId!=cellCount(); cellId++) {
            const double scaling =
                (normalizationMethod==NormalizationMethod::L1) ?
                    sums[cellId].sum1 :
                    sqrt(sums[cellId].sum2);
            if(scaling != 0.) {
                const float factor = float(1./scaling);
                for(GeneId geneId=0; geneId!=geneCount(); geneId++) {
                    v[geneId][cellId] *= factor;
                }
            }
        }
    }

}
