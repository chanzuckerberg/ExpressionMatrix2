#ifndef CZI_EXPRESSION_MATRIX2_EXPRESSION_MATRIX_SUBSET_HPP
#define CZI_EXPRESSION_MATRIX2_EXPRESSION_MATRIX_SUBSET_HPP



/*******************************************************************************

Class ExpressionMatrixSubset is used to store expression counts for a
subset of cells and a subset of genes.

It uses temporary files for all of its mapped files, so it is not persistent.

*******************************************************************************/

#include "CellSets.hpp"
#include "GeneSet.hpp"
#include "Ids.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

#include "utility.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class ExpressionMatrixSubset;
    }
}

class ChanZuckerberg::ExpressionMatrix2::ExpressionMatrixSubset {
public:

    // To create the ExpressionMatrixSubset we need to specify
    // the name of the directory that will contain the temporary files,
    // plus the GeneSet and CellSet to be used
    // and the expression counts for the global expression matrix.
    using CellExpressionCounts = MemoryMapped::VectorOfVectors<pair<GeneId, float>, uint64_t>;
    ExpressionMatrixSubset(
        const string& directoryName,
        const GeneSet& geneSet,
        const CellSet& cellSet,
        const CellExpressionCounts& globalExpressionCounts);

    // The set of genes used by this ExpressionMatrixSubset.
    // Indexed by the local GeneId, contains the global GeneId
    // corresponding to each local GeneId.
    // Stored in order of increasing GeneId, so the mapping
    // from local GeneId to global GeneId is strictly increasing.
    const GeneSet& geneSet;

    // The set of cells used by this ExpressionMatrixSubset.
    // Indexed by the local CellId, contains the global CellId
    // corresponding to each local CellId.
    // Stored in order of increasing CellId, so the mapping
    // from local CellId to global CellId is strictly increasing.
    const CellSet& cellSet;

    // The expression counts for each cell. Stored in sparse format,
    // each with the local GeneId it corresponds to.
    // For each cell, they are stored sorted by increasing GeneId.
    // This is indexed by the local CellId.
    // Local GeneId and CellId are indices in the GeneSet and CellSet in use
    // for this ExpressionbMatrixSubset.
    // This is similar to the cellExpressionCounts of class ExpressionMatrix
    // except that it uses local GeneId's and CellId's instead of global ones.
    CellExpressionCounts cellExpressionCounts;
};

#endif

