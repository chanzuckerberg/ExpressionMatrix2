// This file contains portions of the implementation of class ExpressionMatrix
// that deal with Locality Sensitive Hashing (LSH).
// LSH is used for efficiently finding pairs of similar cells.
// See Chapter 3 of Leskovec, Rajaraman, Ullman, Mining of Massive Datasets,
// Cambridge University Press, 2014, also freely downloadable here:
// http://www.mmds.org/#ver21
// and in particular sections 3.4 through 3.7.

// The similarity between two cells can be written as the cosine of the vector
// between expression counts for the two cells, linearly scaled to zero mean
// and unit variance. As described in section 3.7.2 of the book referenced above,
// we use random unit vectors in gene space to generate LSH functions
// for the cosine similarity.
// These are vectors of dimension equal to the number of genes,
// and with unit L2-norm (the sum of the square if the components is 1).
// These vectors are organized by band and row (see section 3.4.1 of the
// book referenced above). There are lshBandCount bands and lshRowCount
// rows per band, for a total lshBandCount*lshRowCount random vectors.
// Each of these vectors defines an hyperplane orthogonal to it.
// As described in section 3.7.2 of the book referenced above,
// each hyperplane provides a function of a locality-sensitive function.


#include "ExpressionMatrix.hpp"
#include "BitSet.hpp"
#include "ExpressionMatrixSubset.hpp"
#include "heap.hpp"
#include "iterator.hpp"
#include "Lsh.hpp"
#include "multipleSetUnion.hpp"
#include "nextPowerOfTwo.hpp"
#include "orderPairs.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/math/constants/constants.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <cmath>
#include "fstream.hpp"
#include <chrono>
#include <numeric>


// Generate the random unit LSH vectors.
// This uses the Marsaglia method, which consists of generating
// a vector with normally distributed components
// with zero mean and unit variance, then normalizing it, described in
// Marsaglia, G. "Choosing a Point from the Surface of a Sphere." Ann. Math. Stat. 43, 645-646, 1972.
// See http://mathworld.wolfram.com/HyperspherePointPicking.html
void ExpressionMatrix::generateLshVectors(
    GeneId geneCount,
    size_t lshBandCount,
    size_t lshRowCount,
    unsigned int seed,
    vector<vector<vector<double> > >& lshVectors    // Indexed by [band][row][geneId]
    )
{
    // Prepare to generate normally vector distributed components.
    using RandomSource = boost::mt19937;
    using NormalDistribution = boost::normal_distribution<>;
    RandomSource randomSource(seed);
    NormalDistribution normalDistribution;
    boost::variate_generator<RandomSource, NormalDistribution> normalGenerator(randomSource, normalDistribution);

    // Allocate space for the LSH vectors.
    lshVectors.resize(lshBandCount, vector<vector<double> >(lshRowCount, vector<double>(geneCount)));



    // Triple loop over bands, rows, components (genes).
    for(size_t band = 0; band < lshBandCount; band++) {
        for(size_t row = 0; row < lshRowCount; row++) {
            for(GeneId geneId = 0; geneId < geneCount; geneId++) {
                lshVectors[band][row][geneId] = normalGenerator();
            }

            // Normalize the vector for this band and row.
            double sum2 = 0.;
            for(const double& x : lshVectors[band][row]) {
                sum2 += x * x;
            }
            const double factor = 1. / sqrt(sum2);
            for(double& x : lshVectors[band][row]) {
                x *= factor;
            }
        }
    }
}



// Approximate computation of the similarity between two cells using
// Locality Sensitive Hashing (LSH).
// Not to be used for code where performance is important,
// because it recomputes the LSH vector every time.
// Note that this uses all genes.
double ExpressionMatrix::computeApproximateLshCellSimilarity(
    size_t lshBandCount,
    size_t lshRowCount,
    unsigned int seed,
    CellId cellId0,
    CellId cellId1) const
{

    // Generate LSH vectors.
    vector<vector<vector<double> > > lshVectors;
    generateLshVectors(geneCount(), lshBandCount, lshRowCount, seed, lshVectors);

    // Compute scalar products of each LSH vector with the expression counts of each cell.
    vector<double> scalarProducts0, scalarProducts1;
    for(size_t band = 0; band < lshBandCount; band++) {
        for(size_t row = 0; row < lshRowCount; row++) {
            scalarProducts0.push_back(computeExpressionCountScalarProduct(cellId0, lshVectors[band][row]));
            scalarProducts1.push_back(computeExpressionCountScalarProduct(cellId1, lshVectors[band][row]));
        }
    }

    // Count the number of scalar products with the same sign.
    size_t sameSignCount = 0;
    for(size_t i = 0; i < scalarProducts0.size(); i++) {
        const double sign0 = std::copysign(1., scalarProducts0[i]);
        const double sign1 = std::copysign(1., scalarProducts1[i]);
        if(sign0 == sign1) {
            ++sameSignCount;
        }
    }

    // Compute the estimated angle.
    const double angle = boost::math::double_constants::pi
        * (1. - double(sameSignCount) / double(lshBandCount * lshRowCount));
    return std::cos(angle);
}



// Approximate computation of the angle between the expression vectors of two cells
// using Locality Sensitive Hashing (LSH).
// The approximate similarity can be computed as the cosine of this angle.
// This recomputes every time the scalar product of the cell expression vector
// with the LSH vectors.
double ExpressionMatrix::computeApproximateLshCellAngle(
    const vector<vector<vector<double> > >& lshVectors,
    CellId cellId0,
    CellId cellId1) const
{
    // Compute scalar products of each LSH vector with the expression counts of each cell.
    vector<pair<double, double> > scalarProducts;
    for(const auto& bandLshVectors : lshVectors) {
        for(const auto& lshVector : bandLshVectors) {
            scalarProducts.push_back(make_pair(
                computeExpressionCountScalarProduct(cellId0, lshVector),
                computeExpressionCountScalarProduct(cellId1, lshVector)));
        }
    }

    // Count the number of scalar products with the same sign.
    size_t sameSignCount = 0;
    for(size_t i = 0; i < scalarProducts.size(); i++) {
        const double sign0 = std::copysign(1., scalarProducts[i].first);
        const double sign1 = std::copysign(1., scalarProducts[i].second);
        if(sign0 == sign1) {
            ++sameSignCount;
        }
    }

    // Compute the estimated angle.
    const size_t totalCount = scalarProducts.size();
    const double angle = boost::math::double_constants::pi * (1. - double(sameSignCount) / double(totalCount));
    return angle;
}



// Approximate computation of the angle between the expression vectors of two cells
// using Locality Sensitive Hashing (LSH).
double ExpressionMatrix::computeApproximateLshCellAngle(
    const BitSet& signature0,
    const BitSet& signature1,
    double bitCountInverse) const
{

    // Count the number of bits where the two signatures disagree.
    size_t oppositeSignCount = countMismatches(signature0, signature1);

    // Compute the estimated angle.
    const double angle = boost::math::double_constants::pi * double(oppositeSignCount) * bitCountInverse;
    return angle;
}



// Compute the scalar product of an LSH vector with the normalized expression counts of a cell.
double ExpressionMatrix::computeExpressionCountScalarProduct(CellId cellId, const vector<double>& v) const
{
    // Sanity check on the vector being passed in.
    CZI_ASSERT(v.size() == geneCount());

    const Cell& cell = cells[cellId];
    const double mean = cell.sum1 / double(geneCount());
    const double sigmaInverse = 1. / sqrt(cell.sum2 / geneCount() - mean * mean);

    // Initialize the scalar product to what it would be
    // if all counts were zero.
    const double zeroContribution = -mean * sigmaInverse;
    double scalarProduct = 0.;
    for(GeneId geneId = 0; geneId < geneCount(); geneId++) {
        scalarProduct += v[geneId] * zeroContribution;
    }

    // Add the contribution of the non-zero expression counts for this cell.
    for(const auto& p : cellExpressionCounts[cellId]) {
        const GeneId geneId = p.first;
        const float& count = p.second;

        // Accumulate into the scalar product.
        scalarProduct += sigmaInverse * double(count) * v[geneId];
    }
    return scalarProduct;
}



// Write a csv file containing, for every pair of cells,
// the exact similarity and the similarity computed using LSH.
// This uses all genes.
void ExpressionMatrix::writeLshSimilarityComparisonSlow(
    size_t lshBandCount,
    size_t lshRowCount,
    unsigned int seed
    ) const
{
    // Generate LSH vectors.
    vector<vector<vector<double> > > lshVectors;
    generateLshVectors(geneCount(), lshBandCount, lshRowCount, seed, lshVectors);

    // Open the csv file
    ofstream csvOut("LshSimilarityComparison.csv");
    csvOut << "Cell0,Cell1,Exact,LSH,Delta\n";


    // Loop over all pairs.
    for(CellId cellId0=0; cellId0!=cellCount()-1; cellId0++) {
        if((cellId0%1) == 0) {
            cout << timestamp << "Working on cell " << cellId0 << " of " << cells.size() << endl;
        }
        for(CellId cellId1=cellId0+1; cellId1!=cellCount(); cellId1++) {
            const double exactSimilarity = computeCellSimilarity(cellId0, cellId1);
            const double approximateAngle = computeApproximateLshCellAngle(lshVectors, cellId0, cellId1);
            const double approximateSimilarity = std::cos(approximateAngle);
            const double delta = approximateSimilarity - exactSimilarity;
            csvOut << cellId0 << ",";
            csvOut << cellId1 << ",";
            csvOut << exactSimilarity << ",";
            csvOut << approximateSimilarity << ",";
            csvOut << delta << "\n";
            csvOut << flush;
        }
    }
}



// This uses all genes.
void ExpressionMatrix::writeLshSimilarityComparison(
    size_t lshBandCount,
    size_t lshRowCount,
    unsigned int seed
    ) const
{
    // Generate LSH vectors.
    vector<vector<vector<double> > > lshVectors;
    cout << timestamp << "Generating the LSH vectors." << endl;
    generateLshVectors(geneCount(), lshBandCount, lshRowCount, seed, lshVectors);

    // Orthogonalize the LSH vectors.
    // This did not seem to give any benefit, so I turned it off to eliminate the
    // dependency on Lapack.
    // cout << timestamp << "Orthogonalizing the LSH vectors." << endl;
    // orthogonalizeLshVectors(lshVectors, lshBandCount*lshRowCount);

    // Compute the LSH signatures of all cells.
    cout << timestamp << "Computing LSH signatures for all cells." << endl;
    vector<BitSet> signatures;
    computeCellLshSignatures(lshVectors, signatures);

    // Open the csv file
    ofstream csvOut("LshSimilarityComparison.csv");
    csvOut << "Cell0,Cell1,Exact,LSH,Delta\n";


    // Loop over all pairs.
    cout << timestamp << "Computing similarities for all cell pairs." << endl;
    const double bitCountInverse = 1. / double(lshBandCount * lshRowCount);
    for(CellId cellId0=0; cellId0!=cellCount()-1; cellId0++) {
        if((cellId0%100) == 0) {
            cout << timestamp << "Working on cell " << cellId0 << " of " << cells.size() << endl;
        }
        for(CellId cellId1=cellId0+1; cellId1!=cellCount(); cellId1++) {
            const double exactSimilarity = computeCellSimilarity(cellId0, cellId1);
            const double approximateAngle = computeApproximateLshCellAngle(signatures[cellId0], signatures[cellId1], bitCountInverse);
            const double approximateSimilarity = std::cos(approximateAngle);
            const double delta = approximateSimilarity - exactSimilarity;
            csvOut << cellId0 << ",";
            csvOut << cellId1 << ",";
            csvOut << exactSimilarity << ",";
            csvOut << approximateSimilarity << ",";
            csvOut << delta << "\n";
            csvOut << flush;
        }
    }
    cout << timestamp << endl;
}




// Given LSH vectors, compute the LSH signature of all cells.
// The LSH signature of a cell is a bit vector with one bit for each of the LSH vectors.
// Each bit is 1 if the scalar product of the cell expression vector
// (normalized to zero mean and unit variance) with the
// the LSH vector is positive, and 0 otherwise.
void ExpressionMatrix::computeCellLshSignatures(
    const vector<vector<vector<double> > >& lshVectors,
    vector<BitSet>& signatures
    ) const
{
    // Count the total number of lsh vectors, assuming that all the bands have the same
    // number of rows.
    CZI_ASSERT(!lshVectors.empty());
    const size_t lshBandCount = lshVectors.size();
    CZI_ASSERT(!lshVectors.front().empty());
    const size_t lshRowCount = lshVectors.front().size();
    const size_t bitCount = lshBandCount * lshRowCount;

    // Check that all lsh vectors have the same length.
    // This is the number of genes we are using to compute cell signtures.
    const GeneId geneCount = GeneId(lshVectors.front().front().size());
    for(const auto& v: lshVectors) {
        CZI_ASSERT(v.size() == geneCount);
    }
    CZI_ASSERT(geneCount == this->geneCount());

    // Compute the sum of the components of each lsh vector.
    // It is needed below to compute the contribution of the
    // expression counts that are zero.
    vector<double> lshVectorsSums(bitCount);
    size_t index = 0;
    for(size_t band = 0; band < lshBandCount; band++) {
        const auto& bandVectors = lshVectors[band];
        CZI_ASSERT(bandVectors.size() == lshRowCount);  // All bands must have the same number of rows.
        for(size_t row = 0; row < lshRowCount; row++, index++) {
            const auto& rowVector = bandVectors[row];
            CZI_ASSERT(rowVector.size() == geneCount);
            lshVectorsSums[index] = std::accumulate(rowVector.begin(), rowVector.end(), 0.);
        }
    }

    // Initialize the signatures to all zero bits.
    signatures.resize(cellCount(), BitSet(bitCount));

    // Vector to hold the scalar products of the normalized expression vector of a cell
    // with each of the LSH vectors.
    vector<double> scalarProducts(bitCount);

    // Loop over all cells.
    for(CellId cellId = 0; cellId < cellCount(); cellId++) {
        if((cellId > 0) && ((cellId % 100) == 0)) {
            cout << timestamp << "Working on cell " << cellId << endl;
        }

        // Compute the mean and standard deviation for this cell.
        const Cell& cell = cells[cellId];
        const double mean = cell.sum1 / double(geneCount);
        const double sigmaInverse = 1. / sqrt(cell.sum2 / geneCount - mean * mean);

        // Initialize the scalar products to what they would be if all counts were zero.
        const double zeroContribution = -mean * sigmaInverse;
        for(size_t index = 0; index < bitCount; index++) {
            scalarProducts[index] = lshVectorsSums[index] * zeroContribution;
        }

        // Add the contributions of the non-zero expression counts for this cell.
        for(const auto& p : cellExpressionCounts[cellId]) {
            const GeneId geneId = p.first;
            const float& count = p.second;
            const double scaledCount = sigmaInverse * double(count);

            // Accumulate into the scalar products.
            size_t index = 0;
            for(size_t band = 0; band < lshBandCount; band++) {
                const auto& bandVectors = lshVectors[band];
                for(size_t row = 0; row < lshRowCount; row++, index++) {
                    CZI_ASSERT(index < scalarProducts.size());
                    scalarProducts[index] += scaledCount * bandVectors[row][geneId];
                }
            }
        }

        // Set to 1 the signature bits corresponding to positive scalar products.
        auto& cellSignature = signatures[cellId];
        for(size_t index = 0; index < bitCount; index++) {
            if(scalarProducts[index] > 0) {
                cellSignature.set(index);
            }
        }
    }
}



// Same as above, but only for a set of cells given in a vector of cell ids (cell set).
// This still uses all the genes.
void ExpressionMatrix::computeCellLshSignatures(
    const vector<vector<vector<double> > >& lshVectors,
    const MemoryMapped::Vector<CellId>& cellSet,
    vector<BitSet>& signatures
    ) const
{
    // Count the total number of lsh vectors, assuming that all the bands have the same
    // number of rows.
    CZI_ASSERT(!lshVectors.empty());
    const size_t lshBandCount = lshVectors.size();
    CZI_ASSERT(!lshVectors.front().empty());
    const size_t lshRowCount = lshVectors.front().size();
    const size_t bitCount = lshBandCount * lshRowCount;

    // Check that all lsh vectors have the same length.
    // This is the number of genes we are using to compute cell signtures.
    const GeneId geneCount = GeneId(lshVectors.front().front().size());
    for(const auto& x: lshVectors) {
        for(const auto& y: x) {
        CZI_ASSERT(y.size() == geneCount);
        }
    }
    CZI_ASSERT(geneCount == this->geneCount());

    // Compute the sum of the components of each lsh vector.
    // It is needed below to compute the contribution of the
    // expression counts that are zero.
    vector<double> lshVectorsSums(bitCount);
    size_t index = 0;
    for(size_t band = 0; band < lshBandCount; band++) {
        const auto& bandVectors = lshVectors[band];
        CZI_ASSERT(bandVectors.size() == lshRowCount);  // All bands must have the same number of rows.
        for(size_t row = 0; row < lshRowCount; row++, index++) {
            const auto& rowVector = bandVectors[row];
            CZI_ASSERT(rowVector.size() == geneCount);
            lshVectorsSums[index] = std::accumulate(rowVector.begin(), rowVector.end(), 0.);
        }
    }

    // Initialize the signatures to all zero bits.
    signatures.resize(cellSet.size(), BitSet(bitCount));

    // Vector to hold the scalar products of the normalized expression vector of a cell
    // with each of the LSH vectors.
    vector<double> scalarProducts(bitCount);

    // Loop over all cells.
    for(CellId localCellId = 0; localCellId < cellSet.size(); localCellId++) {
        if((localCellId > 0) && ((localCellId % 100) == 0)) {
            cout << timestamp << "Working on cell " << localCellId << " of " << cellSet.size() << endl;
        }

        // Compute the mean and standard deviation for this cell.
        const CellId globalCellId = cellSet[localCellId];
        const Cell& cell = cells[globalCellId];
        const double mean = cell.sum1 / double(geneCount);
        const double sigmaInverse = 1. / sqrt(cell.sum2 / geneCount - mean * mean);

        // Initialize the scalar products to what they would be if all counts were zero.
        const double zeroContribution = -mean * sigmaInverse;
        for(size_t index = 0; index < bitCount; index++) {
            scalarProducts[index] = lshVectorsSums[index] * zeroContribution;
        }

        // Add the contributions of the non-zero expression counts for this cell.
        for(const auto& p : cellExpressionCounts[globalCellId]) {
            const GeneId geneId = p.first;
            const float& count = p.second;
            const double scaledCount = sigmaInverse * double(count);

            // Accumulate into the scalar products.
            size_t index = 0;
            for(size_t band = 0; band < lshBandCount; band++) {
                const auto& bandVectors = lshVectors[band];
                for(size_t row = 0; row < lshRowCount; row++, index++) {
                    CZI_ASSERT(index < scalarProducts.size());
                    scalarProducts[index] += scaledCount * bandVectors[row][geneId];
                }
            }
        }

        // Set to 1 the signature bits corresponding to positive scalar products.
        auto& cellSignature = signatures[localCellId];
        for(size_t index = 0; index < bitCount; index++) {
            if(scalarProducts[index] > 0) {
                cellSignature.set(index);
            }
        }
    }
}



// Same as above, but using a subset of gene and cells.
void ExpressionMatrix::computeCellLshSignatures(
    const ExpressionMatrixSubset& expressionMatrixSubset,
    const vector<vector<vector<double> > >& lshVectors,
    vector<BitSet>& signatures
    )
{
    // Sanity check on the LSH vectors.
    CZI_ASSERT(!lshVectors.empty());
    const size_t lshBandCount = lshVectors.size();          // The number of LSH "bands".
    CZI_ASSERT(!lshVectors.front().empty());
    const size_t lshRowCount = lshVectors.front().size();   // The number of LSH "rows" in each band.
    const GeneId geneCount = expressionMatrixSubset.geneCount();
    for(const auto& band: lshVectors) {
        CZI_ASSERT(band.size() == lshRowCount);
        for(const auto& lshVector: band) {
            CZI_ASSERT(lshVector.size() == geneCount);
        }
    }

    // The number of bits in each signature vector equals the
    // total number of LSH vectors, that is, the total number
    // of rows in all bands.
    const size_t bitCount = lshBandCount * lshRowCount;

    // Compute the sum of the components of each lsh vector.
    // It is needed below to compute the contribution of the
    // expression counts that are zero.
    vector<double> lshVectorsSums(bitCount);
    size_t index = 0;
    for(size_t band = 0; band < lshBandCount; band++) {
        const auto& bandVectors = lshVectors[band];
        CZI_ASSERT(bandVectors.size() == lshRowCount);  // All bands must have the same number of rows.
        for(size_t row = 0; row < lshRowCount; row++, index++) {
            const auto& rowVector = bandVectors[row];
            CZI_ASSERT(rowVector.size() == geneCount);
            lshVectorsSums[index] = std::accumulate(rowVector.begin(), rowVector.end(), 0.);
        }
    }

    // Initialize the signatures to all zero bits.
    signatures.resize(expressionMatrixSubset.cellCount(), BitSet(bitCount));

    // Vector to hold the scalar products of the normalized expression vector of a cell
    // with each of the LSH vectors.
    vector<double> scalarProducts(bitCount);



    // Loop over all cells in our expressionmatrix subset..
    for(CellId localCellId = 0; localCellId < expressionMatrixSubset.cellCount(); localCellId++) {
        if((localCellId > 0) && ((localCellId % 100) == 0)) {
            cout << timestamp << "Working on cell " << localCellId;
            cout << " of " << expressionMatrixSubset.cellCount() << endl;
        }

        // The global CellId should not be needed because both the ExpressionMatrixSubset
        // and the signatures object are indexed by the localCellId.
        // const CellId globalCellId = expressionMatrix.cellSet[localCellId];

        // Compute the mean and standard deviation for this cell.
        const ExpressionMatrixSubset::Sum& sum = expressionMatrixSubset.sums[localCellId];
        const double mean = sum.sum1 / double(geneCount);
        const double sigmaInverse = 1. / sqrt(sum.sum2 / geneCount - mean * mean);

        // Initialize the scalar products to what they would be if all counts were zero.
        const double zeroContribution = -mean * sigmaInverse;
        for(size_t index = 0; index < bitCount; index++) {
            scalarProducts[index] = lshVectorsSums[index] * zeroContribution;
        }

        // Add the contributions of the non-zero expression counts for this cell.
        for(const auto& p : expressionMatrixSubset.cellExpressionCounts[localCellId]) {
            const GeneId localGeneId = p.first;
            const float& count = p.second;
            const double scaledCount = sigmaInverse * double(count);

            // Accumulate into the scalar products.
            size_t index = 0;
            for(size_t band = 0; band < lshBandCount; band++) {
                const auto& bandVectors = lshVectors[band];
                for(size_t row = 0; row < lshRowCount; row++, index++) {
                    CZI_ASSERT(index < scalarProducts.size());
                    scalarProducts[index] += scaledCount * bandVectors[row][localGeneId];
                }
            }
        }

        // Set to 1 the signature bits corresponding to positive scalar products.
        auto& cellSignature = signatures[localCellId];
        for(size_t index = 0; index < bitCount; index++) {
            if(scalarProducts[index] > 0) {
                cellSignature.set(index);
            }
        }
    }
}



// Write to a csv file statistics of the cell LSH signatures..
void ExpressionMatrix::writeLshSignatureStatistics(size_t bitCount, const vector<BitSet>& signatures) const
{
    ofstream csvOut("LshSignatureStatistics.csv");
    csvOut << "Bit,Set,Unset,Total\n";

    for(size_t i = 0; i < bitCount; i++) {

        // Count the number of cells that have this bit set.
        size_t setCount = 0;
        for(CellId cellId = 0; cellId < cellCount(); cellId++) {
            if(signatures[cellId].get(i)) {
                ++setCount;
            }
        }
        const size_t unsetCount = cellCount() - setCount;

        csvOut << i << "," << setCount << "," << unsetCount << "," << cellCount() << "\n";

    }

}




// Find similar cell pairs by looping over all pairs
// and using an LSH approximation to compute the similarity between two cells.
// See the beginning of ExpressionMatrixLsh.cpp for more information.
// Like findSimilarPairs0, this is also O(N**2) slow. However
// the coefficient of the N**2 term is much lower (around 15 ns/pair when nothing gets stored), at a cost of
// additional O(N) work (typically 30 ms per cell for lshCount=1024).
// As a result, this can be much faster for large numbers of cells.
// The error of the approximation is controlled by lshCount.
// The maximum standard deviation of the computed similarity is (pi/2)/sqrt(lshCount),
// or about 0.05 for lshCount=1024.
// The standard deviation decreases as the similarity increases. It becomes
// zero when the similarity is 1. For similarity 0.5, the standard deviation is 82%
// of the standard deviation at similarity 0.
// THIS IS THE OLD VERSION THAT USES ALL THE GENES.
// THE NEW VERSION BELOW ONLY TAKES INTO ACOUNT GENES IN A SPECIFIED GENE SET.
void ExpressionMatrix::findSimilarPairs1Old(
    const string& cellSetName,  // The name of the cell set to be used.
    const string& name,         // The name of the SimilarPairs object to be created.
    size_t k,                   // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold, // The minimum similarity for a pair to be stored.
    size_t lshCount,            // The number of LSH functions (hyperplanes) to be used.
    unsigned int seed           // The seed used to generate the random hyperplanes.
    )
{
    // Sanity check.
    CZI_ASSERT(similarityThreshold <= 1.);

    // For now this uses the AllGenes gene set.
    const string geneSetName = "AllGenes";
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;

    // Locate the cell set.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    if(cellSet.size() == 0) {
        cout << "Cell set " << cellSetName << " is empty. Skipping findSimilarPairs1." << endl;
        return;
    }

    // Generate LSH vectors.
    const size_t lshBandCount = lshCount;
    const size_t lshRowCount = 1;
    vector<vector<vector<double> > > lshVectors;
    generateLshVectors(geneCount(), lshBandCount, lshRowCount, seed, lshVectors);

    // Compute the LSH signatures of all cells.
    cout << timestamp << "Computing LSH signatures for all cells." << endl;
    vector<BitSet> signatures;
    computeCellLshSignatures(lshVectors, cellSet, signatures);

    // Write signature statistics to a csv file.
    // writeLshSignatureStatistics(lshCount, signatures);

    // Create the SimilarPairs object where we will store the pairs.
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + name, k, geneSet, cellSet);

    // Compute the angle threshold (in radians) corresponding to this similarity threshold.
    const double angleThreshold = std::acos(similarityThreshold);

    // Compute the corresponding threshold on the number of mismatching
    // signature bits.
    const size_t mismatchCountThreshold = min(lshCount, size_t(double(lshCount) * angleThreshold / boost::math::double_constants::pi));

    // Compute the similarity (cosine of the angle) corresponding to every number of mismatching bits
    // up to the threshold.
    vector<double> similarityTable(mismatchCountThreshold + 1);
    for(size_t mismatchingBitCount = 0;
        mismatchingBitCount <= mismatchCountThreshold; mismatchingBitCount++) {
        const double angle = double(mismatchingBitCount) *
            boost::math::double_constants::pi / double(lshCount);
        CZI_ASSERT(mismatchingBitCount < similarityTable.size());
        similarityTable[mismatchingBitCount] = std::cos(angle);
    }



    // Loop over all pairs. This is much faster than findSimilarPairs0
    // (around 15 ns per pair when using LSH vectors of 1024 bits), but
    // still scales like the square of the number of cells in the cell set.
    cout << timestamp << "Begin computing similarities for all cell pairs." << endl;
    for(CellId localCellId0=0; localCellId0!=cellSet.size()-1; localCellId0++) {
        if(localCellId0>0 && ((localCellId0%10000) == 0)) {
            cout << timestamp << "Working on cell " << localCellId0 << " of " << cellSet.size() << endl;
        }
        const BitSet& signature0 = signatures[localCellId0];
        for(CellId localCellId1=localCellId0+1; localCellId1!=cellSet.size(); localCellId1++) {
            const BitSet& signature1 = signatures[localCellId1];

            // Count the number of bits where the signatures of these two cells disagree.
            size_t mismatchingBitCount = countMismatches(signature0, signature1);

            // If the similarity is sufficient, pass it to the SimilarPairs container,
            // which will make the decision whether to store it, depending on the
            // number of pairs already stored for cellId0 and cellId1.
            if(mismatchingBitCount <= mismatchCountThreshold) {
                CZI_ASSERT(mismatchingBitCount < similarityTable.size());
                similarPairs.add(localCellId0, localCellId1, similarityTable[mismatchingBitCount]);
            }
        }
    }



    // Sort the similar pairs for each cell by decreasing similarity.
    cout << timestamp << "Sorting pairs." << endl;
    similarPairs.sort();
    cout << timestamp << "Done sorting pairs." << endl;

}



// Find similar cell pairs by looping over all pairs
// and using an LSH approximation to compute the similarity between two cells.
// See the beginning of ExpressionMatrixLsh.cpp for more information.
// Like findSimilarPairs0, this is also O(N**2) slow. However
// the coefficient of the N**2 term is much lower (around 15 ns/pair when nothing gets stored), at a cost of
// additional O(N) work (typically 30 ms per cell for lshCount=1024).
// As a result, this can be much faster for large numbers of cells.
// The error of the approximation is controlled by lshCount.
// The maximum standard deviation of the computed similarity is (pi/2)/sqrt(lshCount),
// or about 0.05 for lshCount=1024.
// The standard deviation decreases as the similarity increases. It becomes
// zero when the similarity is 1. For similarity 0.5, the standard deviation is 82%
// of the standard deviation at similarity 0.
void ExpressionMatrix::findSimilarPairs1(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed               // The seed used to generate the LSH vectors.
    )
{
    // Sanity check.
    CZI_ASSERT(similarityThreshold <= 1.);

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Create the expression matrix subset for this gene set and cell set.
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset-" + similarPairsName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);

    // Generate LSH vectors.
    const size_t lshBandCount = lshCount;
    const size_t lshRowCount = 1;
    vector<vector<vector<double> > > lshVectors;
    generateLshVectors(geneSet.size(), lshBandCount, lshRowCount, seed, lshVectors);

    // Compute the LSH signatures of all cells.
    cout << timestamp << "Computing LSH signatures for all cells." << endl;
    vector<BitSet> signatures;
    computeCellLshSignatures(expressionMatrixSubset, lshVectors, signatures);

    // Compute the angle threshold (in radians) corresponding to this similarity threshold.
    const double angleThreshold = std::acos(similarityThreshold);

    // Compute the threshold on the number of mismatching
    // signature bits.
    size_t mismatchCountThreshold = lshCount;
    if(similarityThreshold > -1.) {
        mismatchCountThreshold = size_t(double(lshCount) * angleThreshold / boost::math::double_constants::pi);
    }

    // Compute the similarity (cosine of the angle) corresponding to every number of mismatching bits
    // up to the threshold.
    vector<double> similarityTable(mismatchCountThreshold + 1);
    for(size_t mismatchingBitCount = 0;
        mismatchingBitCount <= mismatchCountThreshold; mismatchingBitCount++) {
        const double angle = double(mismatchingBitCount) *
            boost::math::double_constants::pi / double(lshCount);
        CZI_ASSERT(mismatchingBitCount < similarityTable.size());
        similarityTable[mismatchingBitCount] = std::cos(angle);
    }

    // Create the SimilarPairs object where we will store the pairs.
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, k, geneSet, cellSet);



    // Loop over all pairs. This is much faster than findSimilarPairs0
    // (around 15 ns per pair when using LSH vectors of 1024 bits), but
    // still scales like the square of the number of cells in the cell set.
    cout << timestamp << "Begin computing similarities for all cell pairs." << endl;
    for(CellId localCellId0=0; localCellId0!=cellCount-1; localCellId0++) {
        if(localCellId0>0 && ((localCellId0%10000) == 0)) {
            cout << timestamp << "Working on cell " << localCellId0 << " of " << cellSet.size() << endl;
        }
        const BitSet& signature0 = signatures[localCellId0];
        for(CellId localCellId1=localCellId0+1; localCellId1!=cellCount; localCellId1++) {
            const BitSet& signature1 = signatures[localCellId1];

            // Count the number of bits where the signatures of these two cells disagree.
            size_t mismatchingBitCount = countMismatches(signature0, signature1);

            // If the similarity is sufficient, pass it to the SimilarPairs container,
            // which will make the decision whether to store it, depending on the
            // number of pairs already stored for cellId0 and cellId1.
            if(mismatchingBitCount <= mismatchCountThreshold) {
                CZI_ASSERT(mismatchingBitCount < similarityTable.size());
                similarPairs.add(localCellId0, localCellId1, similarityTable[mismatchingBitCount]);
            }
        }
    }



    // Sort the similar pairs for each cell by decreasing similarity.
    cout << timestamp << "Sorting pairs." << endl;
    similarPairs.sort();
    cout << timestamp << "Done sorting pairs." << endl;
}



#if 0
// Orthogonalize the LSH vectors in groups of k.
// Ji et al. (2012) have shown that orthogonalization of
// groups of LSH vectors can result in significant reduction of the
// variance of the distance estimate provided by LSH.
// See J. Ji, J. Li, S. Yan, B. Zhang, and Q. Tian,
// Super-Bit Locality-Sensitive Hashing, In NIPS, pages 108–116, 2012.
// https://pdfs.semanticscholar.org/64d8/3ccbcb1d87bfafee57f0c2d49043ee3f565b.pdf
// This did not seem to give any benefit, so I turned it off to eliminate the
// dependency on Lapack.
void ExpressionMatrix::orthogonalizeLshVectors(
    vector< vector< vector<double> > >& lshVectors,
    size_t k
) const
{

    // Create a vector of pair(band, row) for all the LSH vectors,
    // We will then partition that in groups of k.
    const size_t lshBandCount = lshVectors.size();
    CZI_ASSERT(lshBandCount > 0);
    const size_t lshRowCount = lshVectors.front().size();
    CZI_ASSERT(lshRowCount > 0);
    vector< pair<size_t, size_t> > lshVectorsInfo;
    for(size_t band=0; band<lshBandCount; band++) {
        CZI_ASSERT(lshVectors[band].size() == lshRowCount);
        for(size_t row=0; row<lshRowCount; row++) {
            CZI_ASSERT(lshVectors[band][row].size() == geneCount());
            lshVectorsInfo.push_back(make_pair(band, row));
        }
    }

    // Matrix to hold a group of k LSH vectors to be orthogonalized.
    // Each vector is a column of the matrix, of length geneCount(),
    // stored contiguously in memory. This is the format required
    // by Lapack routine dgeqrf.
    vector<double> a(geneCount() * k);

    // Work vectors for Lapack calls.
    vector<double> tau(k);
    const size_t lWork = k * 128;
    vector<double> work(lWork);

    // Loop over groups of k LSH vectors.
    // Orthogonalize each group.
    for(size_t begin=0; begin<lshVectorsInfo.size(); begin+=k) {
        const size_t end = min(begin+k, lshVectorsInfo.size());

        // The number of vectors we are orthogonalizing in this iteration.
        // Always equal to k, except possibly for the last group.
        const size_t n = end - begin;

        // Copy the LSH vectors into the columns of matrix a.
        size_t aIndex = 0ULL;
        for(size_t i=begin; i!=end; i++) {
            const auto& p = lshVectorsInfo[i];
            const size_t band = p.first;
            const size_t row = p.second;
            const vector<double>& v = lshVectors[band][row];
            CZI_ASSERT(v.size() == geneCount());
            copy(v.begin(), v.end(), a.begin()+aIndex);
            aIndex += geneCount();
        }
        CZI_ASSERT(aIndex == n * geneCount());

        // Use Lapack to orthogonalize the vectors stored in a.
        // This uses a Householder QR factorization, which is
        // numerically much more stable than the simple
        // Gram-Schmidt orthogonalization procedure.
        int info = 0;
        dgeqrf_(int(geneCount()), int(n), &a.front(), int(geneCount()), &tau.front(), &work.front(), int(lWork), info);
        CZI_ASSERT(info == 0);
        dorgqr_(int(geneCount()), int(n), int(n), &a.front(), int(geneCount()), &tau.front(), &work.front(), int(lWork), info);
        CZI_ASSERT(info == 0);

        // Copy back the columns of matrix a,
        // which now contain the orthogonalized LSH vectors.
        aIndex = 0ULL;
        for(size_t i=begin; i!=end; i++) {
            const auto& p = lshVectorsInfo[i];
            const size_t band = p.first;
            const size_t row = p.second;
            vector<double>& v = lshVectors[band][row];
            CZI_ASSERT(v.size() == geneCount());
            copy(a.begin()+aIndex, a.begin()+aIndex+geneCount(), v.begin());
            aIndex += geneCount();
        }
        CZI_ASSERT(aIndex == n * geneCount());
    }
}
#endif



// Analyze the quality of a set of similar pairs.
void ExpressionMatrix::analyzeSimilarPairs(
    const string& similarPairsName,
    double csvDownsample) const
{
    // Open the SimilarPairs object we want to analyze.
    const SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, true);
    const GeneSet& geneSet = similarPairs.getGeneSet();
    const CellSet& cellSet = similarPairs.getCellSet();
    const CellId cellCount = CellId(cellSet.size());

    // Create the expression matrix subset to be used
    // for exact similarity computations.
    const string expressionMatrixSubsetName = directoryName + "/tmp-ExpressionMatrixSubset-" + similarPairsName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);

    // Open the output csv file.
    ofstream csvOut(similarPairsName + "-analysis.csv");
    csvOut << "GlobalCellId0,GlobalCellId1,ExactSimilarity,StoredSimilarity\n";

    // Statistics for bins of exact similarity values.
    const size_t binCount = 200;
    const double binWidth = 2. / binCount;
    vector<size_t> sum0(binCount, 0);
    vector<double> sum1(binCount, 0.);
    vector<double> sum2(binCount, 0.);


    // Random number generator used for downsampling
    using RandomSource = boost::mt19937;
    using UniformDistribution = boost::uniform_01<>;
    const int seed = 231;
    RandomSource randomSource(seed);
    UniformDistribution uniformDistribution;
    boost::variate_generator<RandomSource, UniformDistribution>
        uniformGenerator(randomSource, uniformDistribution);

    // Loop over all values stored in the SimilarPairs object.
    size_t pairCount = 0;
    size_t csvPairCount = 0;
    for(CellId localCellId0=0; localCellId0<cellCount; localCellId0++) {
        if((localCellId0 % 100)==0) {
            cout << timestamp << "Working on cell " << localCellId0 << " of " << cellCount << endl;
        }
        const CellId globalCellId0 = cellSet[localCellId0];
        for(const auto& p: similarPairs[localCellId0]) {
            ++pairCount;
            const CellId localCellId1 = p.first;
            const CellId globalCellId1 = cellSet[localCellId1];
            const float& storedSimilarity = p.second;
            const double exactSimilarity = expressionMatrixSubset.
                computeCellSimilarity(localCellId0, localCellId1);

            // Update statistics.
            const double delta = storedSimilarity - exactSimilarity;
            const size_t bin = size_t(floor((exactSimilarity+1.) / binWidth));
            CZI_ASSERT(bin < binCount);
            ++(sum0[bin]);
            sum1[bin] += delta;
            sum2[bin] += delta*delta;

            // Write to csv output (with downsampling).
            if(uniformGenerator() < csvDownsample) {
                ++csvPairCount;
                csvOut << globalCellId0 << ",";
                csvOut << globalCellId1 << ",";
                csvOut << exactSimilarity << ",";
                csvOut << storedSimilarity << "\n";
            }
        }
    }

   cout << "Total number of ordered cell pairs: " << size_t(cellCount) * size_t(cellCount-1) << endl;
   cout << "Number of cell pairs with a stored similarity value: "<< pairCount << endl;
   cout << "Number of cell pairs written to csv file: " << csvPairCount << endl;


   // Compute average and standard deviation of the error for each bin.
   ofstream statsOut(similarPairsName + "-analysis-statistics.csv");
   statsOut << "Similarity,Bias,Rms\n";
   for(size_t bin=0; bin<binCount; bin++) {
       if(sum0[bin] < 2) {
           continue;
       }
       const double similarity = (double(bin) + 0.5) * binWidth - 1.;
       const double s0 = double(sum0[bin]);
       const double s1 = sum1[bin];
       const double s2 = sum2[bin];
       const double average = s1 / s0;
       const double sigma = sqrt(s2 / s0); // Sigma around 0.
       statsOut << similarity << ",";
       statsOut << average << ",";
       statsOut << sigma << "\n";
   }

}



// Find similar cell pairs by looping over all pairs
// and using an LSH approximation to compute the similarity between two cells.
// This is a newer replacement for findSimilarPairs3.
// It is written using class Lsh.
void ExpressionMatrix::findSimilarPairs3(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed               // The seed used to generate the LSH vectors.
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs3 begins." << endl;

    // Sanity check.
    CZI_ASSERT(similarityThreshold <= 1.);

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Create the expression matrix subset for this gene set and cell set.
    cout << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset-" + similarPairsName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);

    // Create the SimilarPairs object where we will store the pairs.
    cout << timestamp << "Initializing SimilarPairs object." << endl;
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, k, geneSet, cellSet);

    // Create the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/tmp-Lsh", expressionMatrixSubset, lshCount, seed);



    // Loop over all pairs. This is much faster than findSimilarPairs0
    // (around 15 ns per pair when using LSH vectors of 1024 bits), but
    // still scales like the square of the number of cells in the cell set.
    // This loops in blocks of cells for better memory locality than a simple
    // loop over cell pairs.
    cout << timestamp << "Begin computing similarities for all cell pairs." << endl;
    const auto t0 = std::chrono::steady_clock::now();
    const CellId blockSize = 64;
    size_t pairCount = 0;
    for(CellId begin0=0; begin0<cellCount; begin0+=blockSize) {
        const CellId end0 = min(begin0+blockSize, cellCount);
        for(CellId begin1=0; begin1<=begin0; begin1+=blockSize) {
            const CellId end1 = min(begin1+blockSize, end0);
            for(CellId cell0=begin0; cell0!=end0; ++cell0) {
                for(CellId cell1=begin1; cell1!=end1 && cell1<cell0; ++cell1) {
                    ++pairCount;

                    // Compute the LSH similarity between these two cells.
                    const double similarity = lsh.computeCellSimilarity(cell0, cell1);

                    // If the similarity is sufficient, pass it to the SimilarPairs container,
                    // which will make the decision whether to store it, depending on the
                    // number of pairs already stored for cellId0 and cellId1.
                    if(similarity > similarityThreshold) {
                        similarPairs.addNoDuplicateCheck(cell0, cell1, similarity);
                    }
                }
            }
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    CZI_ASSERT(pairCount == size_t(cellCount)*(size_t(cellCount-1))/2);
    cout << "Time for all pairs: " << t01 << " s." << endl;
    cout << "Time per pair: " << t01/(0.5*double(cellCount)*double(cellCount-1)) << " s." << endl;


    // Sort the similar pairs for each cell by decreasing similarity.
    cout << timestamp << "Sorting similar pairs." << endl;
    similarPairs.sort();
    cout << timestamp << "ExpressionMatrix::findSimilarPairs3 ends." << endl;

    lsh.remove();

}



// Same as findSimilarPairs3, but without storing anything.
// Used for benchmarking.
void ExpressionMatrix::findSimilarPairs3Benchmark(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed               // The seed used to generate the LSH vectors.
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs3Benchmark begins." << endl;

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Create the expression matrix subset for this gene set and cell set.
    cout << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset";
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);

    // Create the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/tmp-Lsh", expressionMatrixSubset, lshCount, seed);



    const auto t0 = std::chrono::steady_clock::now();
    const CellId blockSize = 64;
    size_t pairCount = 0;
    double sum = 0.;
    for(CellId begin0=0; begin0<cellCount; begin0+=blockSize) {
        const CellId end0 = min(begin0+blockSize, cellCount);
        for(CellId begin1=0; begin1<=begin0; begin1+=blockSize) {
            const CellId end1 = min(begin1+blockSize, end0);
            for(CellId cell0=begin0; cell0!=end0; ++cell0) {
                for(CellId cell1=begin1; cell1!=end1 && cell1<cell0; ++cell1) {
                    ++pairCount;
                    const double similarity = lsh.computeCellSimilarity(cell0, cell1);
                    sum += similarity;
                }
            }
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    CZI_ASSERT(pairCount == size_t(cellCount)*(size_t(cellCount-1))/2);


    cout << "Time for all pairs: " << t01 << " s." << endl;
    cout << "Time per pair: " << t01/(0.5*double(cellCount)*double(cellCount-1)) << " s." << endl;
    cout << "Average similarity: " << sum / (0.5*double(cellCount)*double(cellCount-1)) << endl;

    lsh.remove();

}


// Used to test improvements to findSimilarPairs3.
void ExpressionMatrix::findSimilarPairs4(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed               // The seed used to generate the LSH vectors.
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs4 begins." << endl;

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Create the expression matrix subset for this gene set and cell set.
    cout << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset-" + similarPairsName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);

    // Create the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/tmp-Lsh", expressionMatrixSubset, lshCount, seed);

    // Temporary storage of pairs for each cell.
    vector< vector< pair<CellId, float> > > tmp(cellCount);
    const size_t tmpStore = 2 * k;
    for(auto& v: tmp) {
        v.reserve(tmpStore);
    }

    // Current similarity threshold for each cell.
    vector<float> cellThreshold(cellCount, float(similarityThreshold));



    // Loop over all pairs. This is much faster than findSimilarPairs0
    // (around 15 ns per pair when using LSH vectors of 1024 bits), but
    // still scales like the square of the number of cells in the cell set.
    // This loops in blocks of cells for better memory locality than a simple
    // loop over cell pairs.
    cout << timestamp << "Begin computing similarities for all cell pairs." << endl;
    const auto t0 = std::chrono::steady_clock::now();
    const CellId blockSize = 64;
    size_t pairCount = 0;
    size_t totalPairCount = size_t(cellCount)*(size_t(cellCount-1))/2;
    size_t blockCount = 0;
    for(CellId begin0=0; begin0<cellCount; begin0+=blockSize) {
        const CellId end0 = min(begin0+blockSize, cellCount);
        for(CellId begin1=0; begin1<=begin0; begin1+=blockSize) {
            if(blockCount>0 && ((blockCount%1000000)==0)) {
                cout << timestamp << "Pair computation ";
                cout << 100.*double(pairCount)/double(totalPairCount);
                cout << "% complete." << endl;
            }
            ++blockCount;
            const CellId end1 = min(begin1+blockSize, end0);
            for(CellId cell0=begin0; cell0!=end0; ++cell0) {
                auto& tmp0 = tmp[cell0];
                for(CellId cell1=begin1; cell1!=end1 && cell1<cell0; ++cell1) {
                    auto& tmp1 = tmp[cell1];
                    ++pairCount;

                    // Compute the LSH similarity between these two cells.
                    const double similarity = lsh.computeCellSimilarity(cell0, cell1);

                    // If the similarity is sufficient, pass it to the SimilarPairs container,
                    // which will make the decision whether to store it, depending on the
                    // number of pairs already stored for cellId0 and cellId1.
                    if(similarity > similarityThreshold) {
                        if(similarity > cellThreshold[cell0]) {
                            tmp0.push_back(make_pair(cell1, similarity));
                            if(tmp0.size() == tmpStore) {
                                keepBest(tmp0, k, OrderPairsBySecondGreater< pair<CellId, float> >());
                                cellThreshold[cell0] = tmp0.back().second;
                            }
                        }
                        if(similarity > cellThreshold[cell1]) {
                            tmp1.push_back(make_pair(cell0, similarity));
                            if(tmp1.size() == tmpStore) {
                                keepBest(tmp1, k, OrderPairsBySecondGreater< pair<CellId, float> >());
                                cellThreshold[cell1] = tmp1.back().second;
                            }
                        }
                    }
                }
            }
        }
    }
    // Keep at most k of each.
    for(auto& tmp0: tmp) {
        if(tmp0.size() > k) {
            keepBest(tmp0, k, OrderPairsBySecondGreater< pair<CellId, float> >());
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    CZI_ASSERT(pairCount == totalPairCount);
    cout << "Time for all pairs: " << t01 << " s." << endl;
    cout << "Time per pair: " << t01/(0.5*double(cellCount)*double(cellCount-1)) << " s." << endl;

    // Store the pairs in a SimilarPairs object.
    cout << timestamp << "Initializing SimilarPairs object." << endl;
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, k, geneSet, cellSet);
    cout << timestamp << "Copying similar pairs." << endl;
    similarPairs.copy(tmp);


    // Sort the similar pairs for each cell by decreasing similarity.
    cout << timestamp << "Sorting similar pairs." << endl;
    similarPairs.sort();
    cout << timestamp << "ExpressionMatrix::findSimilarPairs4 ends." << endl;

    lsh.remove();

}



// Replacement for findSimilarPairs2.
// Find similar cell pairs using LSH, without looping over all pairs.
// This uses slices of each LSH signature.
// The number of bits in each slice is lshSliceLength.
// If lshSliceLength is 0, it is set to the base 2 log of the number of cells
// (rounded up).
// For each slice we store the cells with each possible value of the signature slice.
void ExpressionMatrix::findSimilarPairs5(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& lshName,          // The name of the Lsh object to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    size_t lshSliceLength,          // The number of bits in each LSH signature slice, or 0 for automatic selection.
    size_t bucketOverflow           // If not zero, ignore buckets larger than this.
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs5 begins." << endl;
    const auto t0 = std::chrono::steady_clock::now();

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Access the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/Lsh-" + lshName);
    if(lsh.cellCount() != cellSet.size()) {
        throw runtime_error("LSH object " + lshName + " has a number of cells inconsistent with cell set " + cellSetName);
    }


    // Find the signature bits corresponding to each slice.
    const size_t sliceCount = lsh.lshCount() / lshSliceLength;
    vector< vector<size_t> > allSlicesBits(sliceCount);
    for(size_t sliceId=0; sliceId<sliceCount; sliceId++) {
        vector<size_t>& sliceBits = allSlicesBits[sliceId];
        const size_t sliceBegin = sliceId * lshSliceLength;
        const size_t sliceEnd = sliceBegin + lshSliceLength;
        for(size_t bit=sliceBegin; bit!=sliceEnd; bit++) {
            sliceBits.push_back(bit);
        }
    }


    // Table to contain, for each signature slice,
    // the cells with each of possible value of the slice.
    // Indexed by [sliceId][sliceValue].
    // All cell ids are local to the cell set we are using.
    vector< vector < vector<CellId> > > tables(sliceCount);



    // Loop over the signature slices.
    // Each generates a new table of cells.
    for(size_t sliceId=0; sliceId<sliceCount; sliceId++) {
        cout << timestamp << "Computing cell table for signature slice " << sliceId << " of " << sliceCount << endl;
        vector < vector<CellId> >& table = tables[sliceId];
        table.resize(1ULL << lshSliceLength);
        const vector<size_t>& sliceBits = allSlicesBits[sliceId];

        // Store each cell id based on the value of this signature slice.
        for(CellId cellId=0; cellId<cellCount; cellId++) {
            const uint64_t signatureSlice = lsh.getSignature(cellId).getBits(sliceBits);
            CZI_ASSERT(signatureSlice < table.size());
            table[signatureSlice].push_back(cellId);
        }
    }
    cout << timestamp << "Computation of cell table completed." << endl;

    // Temporary storage of pairs for each cell.
    vector< vector< pair<CellId, float> > > tmp(cellCount);


    // Loop over cells.
    // For each cell, compute the union of all the table vectors
    // this cell belongs to. This is the set of candidate neighbors
    // for this cell.
    size_t fullCellCount = 0;
    vector<CellId> candidates;
    vector< const vector<CellId>* > setsToUnion;
    vector< pair<CellId, float> > cellNeighbors;    // The neighbors of a single cell.
    size_t totalCandidateCount = 0;
    for(CellId cellId0=0; cellId0<cellCount; cellId0++) {
        if((cellId0 % 10000)==0) {
            cout << timestamp << "Find neighbors for cell " << cellId0 << " begins." << endl;
        }
        // const auto t0 = std::chrono::steady_clock::now();

        // Find the candidates.
        // This can be made faster using a heap
        // to compute the union.
        candidates.clear();
        setsToUnion.clear();
        // size_t totalCountToUnion = 0;
        for(size_t sliceId=0; sliceId<sliceCount; sliceId++) {
            const uint64_t signatureSlice = lsh.getSignature(cellId0).getBits(allSlicesBits[sliceId]);
            const auto& bucket = tables[sliceId][signatureSlice];
            if(bucketOverflow==0 || bucket.size()<=bucketOverflow) {
                setsToUnion.push_back(&tables[sliceId][signatureSlice]);
                // totalCountToUnion += tables[sliceId][signatureSlice].size();
            }
#if 0
            cout << "Bucket for cell " << cellId0 << " slice " << sliceId << ": ";
            copy(tables[sliceId][signatureSlice].begin(), tables[sliceId][signatureSlice].end(),
                ostream_iterator<CellId>(cout, " "));
            cout << endl;
#endif
        }
        multipleSetUnion(setsToUnion, candidates);
        totalCandidateCount += candidates.size();
        // const auto t1 = std::chrono::steady_clock::now();

        // Check each of the candidates.
        cellNeighbors.clear();
        for(const CellId cellId1: candidates) {
            if(cellId1 == cellId0) {
                continue;
            }
            const double similarity = lsh.computeCellSimilarity(cellId0, cellId1);
            if(similarity > similarityThreshold) {
                cellNeighbors.push_back(make_pair(cellId1, float(similarity)));
            }
        }
        // const auto t2 = std::chrono::steady_clock::now();

#if 0
        cout << "cellNeighbors before keepBest " << cellId0 << endl;
        for(const auto& p: cellNeighbors) {
            cout << p.first << " " << p.second << "\n";
        }
#endif

        // Store the pairs we found, keeping only the k best.
        // const size_t goodCandidatesCount = cellNeighbors.size();
        keepBest(cellNeighbors, k, OrderPairsBySecondGreater< pair<CellId, float> >());
        // const auto t3 = std::chrono::steady_clock::now();
#if 0
        cout << "cellNeighbors after keepBest" << cellId0  << endl;
        for(const auto& p: cellNeighbors) {
            cout << p.first << " " << p.second << "\n";
        }
#endif
        tmp[cellId0] = cellNeighbors;
        // const auto t4 = std::chrono::steady_clock::now();

        if(cellNeighbors.size() == k) {
            ++fullCellCount;
        }

#if 0
        if(cellId0 > 10000) {
            cout << cellId0 << " ";
            cout << 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count()) << " ";
            cout << 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count()) << " ";
            cout << 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count()) << " ";
            cout << 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count()) << " ";
            cout << setsToUnion.size() << " " << totalCountToUnion << " ";
            cout << candidates.size() << " " << goodCandidatesCount << " " << cellNeighbors.size() << "\n";
        }
#endif
    }
    cout << "Average number of candidates per cell is " << double(totalCandidateCount)/cellCount << endl;
    cout << "Number of cells with " << k << " neighbors  is " << fullCellCount << endl;

    // Store the pairs in a SimilarPairs object.
    cout << timestamp << "Initializing SimilarPairs object." << endl;
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, k, geneSet, cellSet);
    cout << timestamp << "Copying similar pairs." << endl;
    similarPairs.copy(tmp);


    // Sort the similar pairs for each cell by decreasing similarity.
    cout << timestamp << "Sorting similar pairs." << endl;
    similarPairs.sort();
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << timestamp << "ExpressionMatrix::findSimilarPairs5 ends. Took " << t01 << " s." << endl;

}



// Find similar cell pairs using LSH and the Charikar algorithm.
// See M. Charikar, "Similarity Estimation Techniques from Rounding Algorithms", 2002,
// section "5. Approximate Nearest neighbor Search in Hamming Space.".
// The Charikar algorithm is for approximate nearest neighbor, but with appropriate
// choices of the algorithm parameters permutationCount and searchCount
// can be used for approximate k nearest neighbors.
// In the Charikar paper, permutationCount is N and searchCount is 2N.
void ExpressionMatrix::findSimilarPairs6(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& lshName,          // The name of the Lsh object to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    size_t permutationCount,        // The number of bit permutations for the Charikar algorithm.
    size_t searchCount,             // The number of cells checked for each cell, in the Charikar algorithm.
    int seed                        // The seed used to randomly generate the bit permutations.
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs6 begins." << endl;
    const bool debug = true;
    const auto t0 = std::chrono::steady_clock::now();

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Access the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/Lsh-" + lshName);
    if(lsh.cellCount() != cellSet.size()) {
        throw runtime_error("LSH object " + lshName + " has a number of cells inconsistent with cell set " + cellSetName);
    }
    const size_t lshCount = lsh.lshCount();


    // Write out the signatures.
    if(debug) {
        cout << "Cell signatures:\n";
        for(size_t i=0; i<lshCount; i++) {
            cout << (i%10);
        }
        cout << "\n";
        for(CellId cellId=0; cellId<cellCount; cellId++) {
            cout << lsh.getSignature(cellId).getString(lshCount) << " " << cellId << "\n";
        }
    }



    // Create the random number generator that will be used to generate
    // the random permutations of the signature bits.
    std::mt19937 randomGenerator(seed);


    // For each of the permutations, we will store:
    // - The permuted signatures, in sorted order.
    // - The corresponding cell ids, in order consistent with the permuted signatures.
    vector< pair<BitSets, vector<CellId> > > permutationData(permutationCount, make_pair(
        BitSets(lsh.cellCount(), lsh.wordCount()),
        vector<CellId>(cellCount)
        ));



    // For each of the permutations, compute permuted/sorted signatures.
    for(size_t permutationId=0; permutationId<permutationCount; permutationId++) {

        // Generate a random permutation of the signature bits.
        vector<uint64_t> bitPermutation(lshCount);
        std::iota(bitPermutation.begin(), bitPermutation.end(), 0ULL);
        std::shuffle(bitPermutation.begin(), bitPermutation.end(), randomGenerator);

        if(debug) {
            cout << "Creating permutation data for permutation " << permutationId << ".\n";
            cout << "Bit permutation:\n";
            for(size_t i=0; i<lshCount; i++) {
                cout << i << " " << bitPermutation[i] << "\n";
            }
        }

        // Compute the permuted signatures for this permutation.
        BitSets permutedSignatures(cellCount, lsh.wordCount());
        for(CellId cellId=0; cellId<cellCount; cellId++) {
            BitSetPointer signature = lsh.getSignature(cellId);
            BitSetPointer permutedSignature = permutedSignatures[cellId];
            permutedSignature.fillUsingPermutation(bitPermutation, signature);
        }

        // Write the permuted signatures.
        if(debug) {
            for(CellId cellId=0; cellId<cellCount; cellId++) {
                cout << permutedSignatures[cellId].getString(lshCount) << " " << cellId << "\n";
            }
        }

        // Sort the permuted signatures lexicographically,
        // Keeping track of the cell ids as they get reordered.
        vector< pair<BitSetPointer, CellId> > table(cellCount);
        for(CellId cellId=0; cellId<cellCount; cellId++) {
            pair<BitSetPointer, CellId>& p = table[cellId];
            p.first = permutedSignatures[cellId];
            p.second = cellId;
        }
        if(debug) {
            cout << "Table before sorting:" << endl;
            for(CellId cellId=0; cellId<cellCount; cellId++) {
                pair<BitSetPointer, CellId>& p = table[cellId];
                cout << p.first.getString(lshCount) << " " << p.second << "\n";
            }
        }
        sort(table.begin(), table.end());
        if(debug) {
            cout << "Table after sorting:" << endl;
            for(CellId i=0; i<cellCount; i++) {
                pair<BitSetPointer, CellId>& p = table[i];
                cout << p.first.getString(lshCount) << " " << p.second << "\n";
            }
        }

        // Store the sorted signatures and corresponding cell ids for this permutation.
        BitSets& thisPermutationBitSets = permutationData[permutationId].first;
        vector<CellId>& thisPermutationCellIds = permutationData[permutationId].second;
        CZI_ASSERT(thisPermutationBitSets.bitSetCount == cellCount);
        CZI_ASSERT(thisPermutationBitSets.wordCount == lsh.wordCount());
        CZI_ASSERT(thisPermutationCellIds.size() == cellCount);
        for(CellId i=0; i<cellCount; i++) {
            cout << "***A " << i << endl;
            pair<BitSetPointer, CellId>& p = table[i];
            cout << "***B " << i << endl;
            thisPermutationBitSets.set(i, p.first);
            cout << "***C " << i << endl;
            thisPermutationCellIds[i] = p.second;
            cout << "***D " << i << endl;
        }

        if(debug) {
            cout << "Permutation data for permutation " << permutationId << ":" << endl;
            for(CellId i=0; i<cellCount; i++) {
                cout << thisPermutationBitSets[i].getString(lshCount) << " " << thisPermutationCellIds[i] << "\n";
            }
        }
    }



    // The rest is not done.
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << timestamp << "ExpressionMatrix::findSimilarPairs6 ends. Took " << t01 << " s." << endl;

    CZI_ASSERT(0);
}



// Compute cell LSH signatures and store them.
void ExpressionMatrix::computeLshSignatures(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& lshName,          // The name of the Lsh object to be created.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed               // The seed used to generate the LSH vectors.
    )
{
    cout << timestamp << "ExpressionMatrix::computeLshSignatures begins." << endl;

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Create the expression matrix subset for this gene set and cell set.
    cout << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset-" + lshName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);

    // Create the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/Lsh-" + lshName, expressionMatrixSubset, lshCount, seed);

    cout << timestamp << "ExpressionMatrix::computeLshSignatures ends." << endl;
}



// Compare two SimilarPairs objects computed using LSH,
// assuming that the first one was computed using a complete
// loop on all pairs (findSimilarPairs4).
void ExpressionMatrix::compareSimilarPairs(
    const string& similarPairsName0,
    const string& similarPairsName1)
{
    // Access the SimilarPairs objects.
    const SimilarPairs similarPairs0(directoryName + "/SimilarPairs-" + similarPairsName0, true);
    const SimilarPairs similarPairs1(directoryName + "/SimilarPairs-" + similarPairsName1, true);

    // Sanity check that the two use the same gene sets and cell sets.
    CZI_ASSERT(similarPairs0.getGeneSet() == similarPairs1.getGeneSet());
    CZI_ASSERT(similarPairs0.getCellSet() == similarPairs1.getCellSet());

    // Loop over cells.
    ofstream csvOut("CompareSimilarPairs.csv");
    csvOut << "CellId,Stored0,Stored1,Lowest0,Lowest1,\n";
    const CellId cellCount = CellId(similarPairs0.getCellSet().size());
    for(CellId cellId=0; cellId<cellCount; cellId++) {
        const auto n0 = similarPairs0.size(cellId);
        const auto n1 = similarPairs1.size(cellId);
        const auto lowest0 = n0 ? ((similarPairs0.end(cellId)-1)->second) : 1.;
        const auto lowest1 = n1 ? ((similarPairs1.end(cellId)-1)->second) : 1.;
        if(n0==n1 && lowest0==lowest1) {
            continue;
        }
        csvOut << cellId << ",";
        csvOut << n0 << ",";
        csvOut << n1 << ",";
        if(n0) {
            csvOut << lowest0;
        }
        csvOut << ",";
        if(n1) {
            csvOut << lowest1;
        }
        csvOut << ",";
        csvOut << "\n";


    }


}


// Analyze the quality of the LSH computation of cell similarity.
void ExpressionMatrix::analyzeLsh(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed,              // The seed used to generate the LSH vectors and to downsample.
    double csvDownsample            // The fraction of pairs that will be included in the output spreadsheet.
    )
{

    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Create the expression matrix subset for this gene set and cell set.
    cout << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset";
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);


    // Create the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/tmp-Lsh", expressionMatrixSubset, lshCount, seed);

    // Random number generator used for downsampling
    using RandomSource = boost::mt19937;
    using UniformDistribution = boost::uniform_01<>;
    RandomSource randomSource(seed);
    UniformDistribution uniformDistribution;
    boost::variate_generator<RandomSource, UniformDistribution>
        uniformGenerator(randomSource, uniformDistribution);

    // Statistics for bins of exact similarity values.
    const size_t binCount = 200;
    const double binWidth = 2. / binCount;
    vector<size_t> sum0(binCount, 0);
    vector<double> sum1(binCount, 0.);
    vector<double> sum2(binCount, 0.);

    // Open the output csv file.
    ofstream csvOut( "Lsh-analysis.csv");
    csvOut << "LocalCellId0,LocalCellId1,GlobalCellId0,GlobalCellId1,ExactSimilarity,LshSimilarity\n";


    // Loop over pairs of cells.
    for(CellId localCellId0=0; localCellId0<cellCount-1; localCellId0++) {
        if((localCellId0%1000) == 0 ) {
            cout << timestamp << "Working on cell " << localCellId0 << " of " << cellCount << endl;
        }
        for(CellId localCellId1=localCellId0+1; localCellId1<cellCount; localCellId1++) {

            // Compute exact similarity for this pair.
            const double exactSimilarity = expressionMatrixSubset.
                computeCellSimilarity(localCellId0, localCellId1);

            // Compute LSH similarity for this pair.
            const double lshSimilarity = lsh.computeCellSimilarity(localCellId0, localCellId1);

            // Update statistics.
            const double delta = lshSimilarity - exactSimilarity;
            const size_t bin = size_t(floor((exactSimilarity+1.) / binWidth));
            CZI_ASSERT(bin < binCount);
            ++(sum0[bin]);
            sum1[bin] += delta;
            sum2[bin] += delta*delta;

            // Write to the output csv file, subject to downsampling.
            if(uniformGenerator() < csvDownsample) {
                csvOut << localCellId0 << ",";
                csvOut << localCellId1 << ",";
                csvOut << cellSet[localCellId0] << ",";
                csvOut << cellSet[localCellId1] << ",";
                csvOut << exactSimilarity << ",";
                csvOut << lshSimilarity << ",\n";
            }
        }

    }



    // Compute average and standard deviation of the error for each bin.
    ofstream statsOut("LSH-analysis-statistics.csv");
    statsOut << "Similarity,Bias,Rms,RmsTheory\n";
    for(size_t bin=0; bin<binCount; bin++) {
        if(sum0[bin] < 2) {
            continue;
        }
        using boost::math::double_constants::pi;
        const double similarity = (double(bin) + 0.5) * binWidth - 1.;
        const double sinTheta = sqrt(1.-similarity*similarity);
        const double theta = std::acos(similarity);
        const double p = 1.- theta / pi;
        const double theoreticalSigma = pi * sinTheta * sqrt(p*(1.-p)/double(lshCount));
        const double s0 = double(sum0[bin]);
        const double s1 = sum1[bin];
        const double s2 = sum2[bin];
        const double average = s1 / s0;
        const double sigma = sqrt(s2 / s0); // Sigma around 0.
        statsOut << similarity << ",";
        statsOut << average << ",";
        statsOut << sigma << ",";
        statsOut << theoreticalSigma << "\n";
    }

    lsh.remove();
}



// Analyze LSH signatures.
void ExpressionMatrix::analyzeLshSignatures(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed              // The seed used to generate the LSH vectors and to downsample.
    )
{
    // Locate the gene set and verify that it is not empty.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = itGeneSet->second;
    if(geneSet.size() == 0) {
        throw runtime_error("Gene set " + geneSetName + " is empty.");
    }

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Create the expression matrix subset for this gene set and cell set.
    cout << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset";
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);


    // Create the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/tmp-Lsh", expressionMatrixSubset, lshCount, seed);


    // Create a map that gives the cells with a given signature.
    map<BitSet, vector<CellId> > signatureMap;
    for(CellId cellId=0; cellId<cellCount; cellId++) {
        const BitSet signature(lsh.getSignature(cellId));
        signatureMap[signature].push_back(cellId);
    }

    // Create a table of signatures ordered by decreasing number of cells.
    vector< pair<BitSet, size_t> > signatureTable;
    for(const auto& p: signatureMap) {
        const BitSet signature = p.first;
        const size_t size = p.second.size();
        signatureTable.push_back(make_pair(signature, size));
    }
#if 0
    sort(signatureTable.begin(), signatureTable.end(),
        OrderPairsBySecondGreater< pair<BitSet, size_t> >());
#endif



    // Write a csv file with one line for each distinct signature.
    {
        ofstream csvOut("Signatures.csv");
        for(const auto& p: signatureTable) {
            const BitSet signature = p.first;
            const size_t size = p.second;
            csvOut << signature.getString(lshCount) << "," << size << "\n";
        }
    }



    vector<size_t> histogram;
    for(const auto& p: signatureMap) {
        const size_t size = p.second.size();
        if(size >= histogram.size()) {
            histogram.resize(size+1, 0);
        }
        ++(histogram[size]);
    }
    ofstream csvOut("Histogram.csv");
    size_t sum = 0;
    for(size_t i=0; i<histogram.size(); i++) {
        const size_t frequency = histogram[i];
        if(frequency) {
            sum += frequency*i;
            csvOut << i << "," << frequency << "," << frequency*i  << "," << sum << "\n";
        }
    }
    cout << flush;


    lsh.writeSignatureStatistics("LshSignatureStatistics.csv");
    lsh.remove();

}
