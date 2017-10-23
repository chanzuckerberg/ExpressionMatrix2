#include "Lsh.hpp"
#include "ExpressionMatrixSubset.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "fstream.hpp"


Lsh::Lsh(
    const ExpressionMatrixSubset& expressionMatrixSubset,
    size_t lshCount,                // Number of LSH hyperplanes
    uint32_t seed                   // Seed to generate LSH hyperplanes.
    )
{
    cout << timestamp << "Generating LSH vectors." << endl;
    generateLshVectors(expressionMatrixSubset.geneCount(), lshCount, seed);

    cout << timestamp << "Computing cell LSH signatures." << endl;
    computeCellLshSignatures(expressionMatrixSubset);
}



// Generate the LSH vectors.
void Lsh::generateLshVectors(
    size_t geneCount,
    size_t lshCount,                // Number of LSH hyperplanes
    uint32_t seed                   // Seed to generate LSH hyperplanes.
)
{
    // Prepare to generate normally vector distributed components.
    using RandomSource = boost::mt19937;
    using NormalDistribution = boost::normal_distribution<>;
    RandomSource randomSource(seed);
    NormalDistribution normalDistribution;
    boost::variate_generator<RandomSource, NormalDistribution> normalGenerator(randomSource, normalDistribution);

    // Allocate space for the LSH vectors.
    lshVectors.resize(geneCount, vector<double>(lshCount, 0.));

    // Sum of the squares of the components of each of the LSH vectors.
    vector<double> normalizationFactor(lshCount, 0.);



    // Loop over genes.
    for(size_t geneId = 0; geneId<lshCount; geneId++) {

        // For this gene, generate the components of all the LSH vectors.
        for(size_t lshVectorId = 0; lshVectorId<lshCount; lshVectorId++) {

            const double x = normalGenerator();
            lshVectors[geneId][lshVectorId] = x;

            // Update the sum of squares for this LSH vector.
            normalizationFactor[lshVectorId] += x*x;
        }
    }

    // Normalize each of the LSH vectors.
    for(auto& f: normalizationFactor) {
        f = 1. / sqrt(f);
    }
    for(size_t geneId = 0; geneId<lshCount; geneId++) {
        for(size_t lshVectorId = 0; lshVectorId<lshCount; lshVectorId++) {
            lshVectors[geneId][lshVectorId] *= normalizationFactor[lshVectorId];
        }
    }

}



// Compute the LSH signatures of all cells in the cell set we are using.
void Lsh::computeCellLshSignatures(const ExpressionMatrixSubset& expressionMatrixSubset)
{
    // Get the number of LSH vectors.
    CZI_ASSERT(!lshVectors.empty());
    const size_t lshCount = lshVectors.front().size();

    // Initialize the cell signatures.
    cout << timestamp << "Initializing cell LSH signatures." << endl;
    const auto cellCount = expressionMatrixSubset.cellCount();
    signatures.resize(cellCount, BitSet(lshCount));

    // Vector to contain, for a single cell, the scalar products of the shifted
    // expression vector for the cell with all of the LSH vectors.
    vector<double> scalarProducts(lshCount);

    // Loop over all the cells in the cell set we are using.
    // The CellId is local to the cell set we are using.
    cout << timestamp << "Computation of cell LSH signatures begins." << endl;
    for(CellId cellId=0; cellId<cellCount; cellId++) {
        if((cellId % 1000) == 0) {
            cout << timestamp << "Working on cell " << cellId << " of " << cellCount << endl;
        }

    }
    cout << timestamp << "Computation of cell LSH signatures ends." << endl;

    CZI_ASSERT(0);
}

