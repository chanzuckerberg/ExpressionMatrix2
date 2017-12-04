// Class to describe an RNA expression matrix.

#ifndef CZI_EXPRESSION_MATRIX2_RAND_INDEX_HPP
#define CZI_EXPRESSION_MATRIX2_RAND_INDEX_HPP



// Computation of the Rand Index and Adjusted Rand Index given a contingency table.
// References:
// - Jorge M. Santos and Mark Embrechts, On the Use of the Adjusted Rand Index
//   as a Metric for Evaluating Supervised Classication,
//   International Conference on Artificial Neural Networks 2009 pp. 175-184,
//   https://pdfs.semanticscholar.org/52d4/8b393f3f838f2370c50af03703eee0bbd669.pdf
// - https://en.wikipedia.org/wiki/Rand_index

#include "algorithm.hpp"
#include "CZI_ASSERT.hpp"
#include "vector.hpp"
#include <numeric>

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        template<class T> void computeRandIndex(
            const vector< vector<T> >& contingencyTable,
            double& randIndex,
            double& adjustedRandIndex)
        {

            // Extract the number of rows and columns of the contingency table.
            const size_t rowCount = contingencyTable.size();
            CZI_ASSERT(rowCount > 0);
            const size_t columnCount = contingencyTable.front().size();
            for(const auto& row: contingencyTable) {
                CZI_ASSERT(row.size() == columnCount);  // All rows must have the same number of columns.
            }

            // Compute row totals.
            vector<T> rowTotals;
            for(const auto& row: contingencyTable) {
                rowTotals.push_back(std::accumulate(row.begin(), row.end(), 0ULL));
            }
            CZI_ASSERT(rowTotals.size() == rowCount);

            // Compute column totals.
            vector<T> columnTotals(columnCount, 0.);
            for(const auto& row: contingencyTable) {
                CZI_ASSERT(row.size() == columnCount);
                for(size_t column=0; column<columnCount; column++) {
                    columnTotals[column] += row[column];
                }
            }
            CZI_ASSERT(columnTotals.size() == columnCount);

            // Compute the total number of elements.
            const T n = std::accumulate(rowTotals.begin(), rowTotals.end(), 0ULL);
            CZI_ASSERT(n == std::accumulate(columnTotals.begin(), columnTotals.end(), 0ULL));
            const double nDouble = double(n);

            // Compute the binomial n over 2.
            const double nBinomial2 = 0.5 * nDouble * (nDouble-1.);

            // Compute a using equation(1) of Santos and Embrechts.
            double a = 0.;
            for(const auto& row: contingencyTable) {
                for(const T& value: row) {
                    const double v = double(value);
                    a += v * (v-1.);
                }
            }
            a /= 2.;

            // Compute b using equation(2) of Santos and Embrechts.
            double b = -a;
            for(const auto rowTotal: rowTotals) {
                const double t = double(rowTotal);
                b += 0.5 * t * (t-1.);
            }

            // Compute c using equation(3) of Santos and Embrechts.
            double c = -a;
            for(const auto columnTotal: columnTotals) {
                const double t = double(columnTotal);
                c += 0.5 * t * (t-1.);
            }

            // Compute d using equation(4) of Santos and Embrechts.
            const double d = nBinomial2 - a - b - c;

            // Compute the Rand Index.
            randIndex = (a + d) / (a + b + c + d);

            // Compute the Adjusted Rand Index.
            const double commonTerm = (a+b)*(a+c) + (c+d)*(b+d);
            const double adjustedRandIndexNumerator = nBinomial2 * (a+d) - commonTerm;
            const double adjustedRandIndexDenominator = nBinomial2*nBinomial2 - commonTerm;
            adjustedRandIndex = adjustedRandIndexNumerator / adjustedRandIndexDenominator;
        }

    }
}


#endif


