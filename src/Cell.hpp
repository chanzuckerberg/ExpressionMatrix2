// Class to describe a single cell corresponding to a row of an RNA expression matrix.

#ifndef CZI_EXPRESSION_MATRIX2_CELL_HPP
#define CZI_EXPRESSION_MATRIX2_CELL_HPP

#include "cstddef.hpp"


namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class Cell;
    }
}


// Class used to store fixed-size information on a single cell.
class ChanZuckerberg::ExpressionMatrix2::Cell {
public:



    // The sum of the expression counts for this cell.
    double sum1;

    // The sum of the squares of the expression counts for this cell.
    double sum2;


    // The L1 norm of the expression counts for this cell.
    // The expression counts cannot be negative.
    // Therefore, the L1 norm equals the sum of the expression counts.
    double norm1() const
    {
        return sum1;
    }

    // The L2 norm of the expression counts for this cell.
    // This equals the square root of the sum of the squares.
    double norm2;

    // Inverses of the L1 and L2 norms.
    double norm1Inverse;
    double norm2Inverse;

    // Quantities used for by computeApproximateCellSimilarity.
    double sum1LargeExpressionCounts;
    double sum2LargeExpressionCounts;


#if 0

    // Additional fields used to compute bounds on approximate cell similarity.
    // None of the bounds turned out to be sufficiently tight, so we are not using these.

    // Quantities used by computeApproximateCellSimilarity.
    // The complete expression counts are stored in ExpressionMatrix::cellExpressionCounts.
    // The large expression counts are stored in ExpressionMatrix::largeCellExpressionCounts.
    // The small expression counts are the one present in cellExpressionCounts but not
    // in largeCellExpressionCounts.
    // Call x and y the complete expression counts for two cells.
    // Call u and v the corresponding large expression counts.
    // The computation of the similarity between two cells requires computing
    // the scalar product x*y.
    // Since all expression counts are non-negative, the scalar product u*v,
    // which is much faster to compute, is a lower bound for x*y.
    // If we define the vector du and dv such that x=u*du, y=v+dv,
    // the maximum error incurred when approximating x*y with u*v can be
    // bounded as follows:
    // x*y = (u*du)*(v+dv) = u*v + u*dv + v*du + du*dv
    // x*y - u*v = u*dv + v*du + du*dv
    // By the triangular inequality:
    // |x*y - u*v| <= |u*dv| + |v*du| + |du*dv|
    // The absolute value of the scalar product of two vectors is <= than the product
    // of the lengths of the vectors (2-norms). Therefore:
    // |x*y - u*v| <= |u|*|dv| + |v|*|du| + |du|*|dv|
    // where |u|, |v|, |du|, and |dv| represents the vector lengths (2-norms).
    // Finally, based on the observation abobve, we know that u*v is a lower bound for x*y.
    // Therefore we can remove the absolute value on the left and write:
    // x*y <= u*v + |u|*|dv| + |v|*|du| + |du|*|dv|
    // This expression gives the upper bound compute by computeApproximateCellSimilarity.
    double largeExpressionCountsNorm2;  // |u| for this cell
    double smallExpressionCountsNorm2;  // |du| for this cell.


    // L2 norms of U and delta U (see "Fast computation of cell similarity",
    // by Paolo Carnevali, dated 5/9/2017.
    double norm2U;
    double norm2DeltaU;

    // Quantities needed to compute the upper bound on page 3 of
    // "Fast computation of cell similarity", by Paolo Carnevali, dated 5/8/2017.
    double largestCount;
    double largestSmallCount;
    double sum1SmallCounts;
#endif

};

#endif

