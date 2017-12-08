#include "Lsh.hpp"
#include "ExpressionMatrixSubset.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/math/constants/constants.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <chrono>
#include "fstream.hpp"



Lsh::Lsh(
    const string& name,             // Name prefix for memory mapped files.
    const ExpressionMatrixSubset& expressionMatrixSubset,
    size_t lshCount,                // Number of LSH hyperplanes
    uint32_t seed                   // Seed to generate LSH hyperplanes.
    )
{
    // Store the Info object.
    info.createNew(name + "-Info");
    info->lshCount = lshCount;
    info->cellCount = expressionMatrixSubset.cellCount();

    // Generate the LSH vectors.
    cout << timestamp << "Generating LSH vectors." << endl;
    generateLshVectors(expressionMatrixSubset.geneCount(), lshCount, seed);

    // Compute cell signatures.
    cout << timestamp << "Computing cell LSH signatures." << endl;
    computeCellLshSignatures(name, expressionMatrixSubset);

    // Compute the similarity table.
    // This is a look up table indexed by the number of mismatching bits.
    // Each entry contgains the similarity corresponding to that number
    // of mismatching bits.
    computeSimilarityTable();
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
    for(size_t geneId = 0; geneId<geneCount; geneId++) {

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
    for(size_t geneId = 0; geneId<geneCount; geneId++) {
        for(size_t lshVectorId = 0; lshVectorId<lshCount; lshVectorId++) {
            lshVectors[geneId][lshVectorId] *= normalizationFactor[lshVectorId];
        }
    }

}



// Compute the LSH signatures of all cells in the cell set we are using.
void Lsh::computeCellLshSignatures(
    const string& name,             // Name prefix for memory mapped files.
    const ExpressionMatrixSubset& expressionMatrixSubset)
{
    // Get the number of LSH vectors.
    CZI_ASSERT(!lshVectors.empty());
    const size_t lshCount = info->lshCount;

    // Compute the number of 64 bit words in each cell signature.
    signatureWordCount = (lshCount-1)/64 + 1;

    // Get the number of genes and cells in the gene set and cell set we are using.
    const auto geneCount = expressionMatrixSubset.geneCount();
    const auto cellCount = expressionMatrixSubset.cellCount();
    CZI_ASSERT(lshVectors.size() == geneCount);

    // Compute the sum of the components of each lsh vector.
    // It is needed below to compute the contribution of the
    // expression counts that are zero.
    vector<double> lshVectorsSums(lshCount, 0.);
    for(GeneId localGeneId=0; localGeneId!=geneCount; localGeneId++) {
        const auto& v = lshVectors[localGeneId]; // Components of all LSH vectors for this gene.
        CZI_ASSERT(v.size() == lshCount);
        for(size_t i=0; i<lshCount; i++) {
            lshVectorsSums[i] += v[i];
        }
    }

    // Initialize the cell signatures.
    cout << timestamp << "Initializing cell LSH signatures." << endl;
    signatures.createNew(name + "-Signatures", cellCount*signatureWordCount);

    // Vector to contain, for a single cell, the scalar products of the shifted
    // expression vector for the cell with all of the LSH vectors.
    vector<double> scalarProducts(lshCount);



    // Loop over all the cells in the cell set we are using.
    // The CellId is local to the cell set we are using.
    cout << timestamp << "Computation of cell LSH signatures begins." << endl;
    const auto t0 = std::chrono::steady_clock::now();
    for(CellId localCellId=0; localCellId<cellCount; localCellId++) {
        if((localCellId % 10000) == 0) {
            cout << timestamp << "Working on cell " << localCellId << " of " << cellCount << endl;
        }

        // Compute the mean of the expression vector for this cell.
        const ExpressionMatrixSubset::Sum& sum = expressionMatrixSubset.sums[localCellId];
        const double mean = sum.sum1 / double(geneCount);

        // If U is one of the LSH vectors, we need to compute the scalar product
        // s = X*U, where X is the cell expression vector, shifted to zero mean:
        // X = x - mean,
        // mean = sum(x)/geneCount (computed above).
        // We get:
        // s = (x-mean)*U = x*U - mean*U = x*U - mean*sum(U)
        // We computed sum(U) above and stored it in lshVectorSums for
        // each of the LSH vectors.
        // Initialize the scalar products for this cell
        // with all of the LSH vectors to -mean*sum(U).
        for(size_t i=0; i<lshCount; i++) {
            scalarProducts[i] = -mean * lshVectorsSums[i];
        }

        // Now add to each scalar product the x*U portion.
        // For performance, the loop over genes is outside,
        // which gives better memory locality.
        // Add the contributions of the non-zero expression counts for this cell.
        for(const auto& p : expressionMatrixSubset.cellExpressionCounts[localCellId]) {
            const GeneId localGeneId = p.first;
            const double count = double(p.second);

            // Add the contribution of this gene to the scalar products.
            const auto& v = lshVectors[localGeneId];
            CZI_ASSERT(v.size() == lshCount);
            for(size_t i=0; i<lshCount; i++) {
                scalarProducts[i] += count * v[i];
            }
        }

        // Set to 1 the signature bits corresponding to positive scalar products.
        BitSetInMemory cellSignature = getSignature(localCellId);
        for(size_t i=0; i<lshCount; i++) {
            if(scalarProducts[i]>0.) {
                cellSignature.set(i);
            }
        }

    }
    const auto t1 = std::chrono::steady_clock::now();
    cout << timestamp << "Computation of cell LSH signatures ends." << endl;
    const size_t nonZeroExpressionCount = expressionMatrixSubset.totalExpressionCounts();
    cout << "Processed " << nonZeroExpressionCount << " non-zero expression counts for ";
    cout << geneCount << " genes and " << cellCount << " cells." << endl;
    cout << "Average number of expression counts per cell  is " << double(nonZeroExpressionCount) / double(cellCount) << endl;
    cout << "Average expression matrix sparsity is " <<
        double(nonZeroExpressionCount) / (double(geneCount) * double(cellCount)) << endl;
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << "Computation of LSH cell signatures took " << t01 << "s." << endl;
    cout << "    Seconds per cell " << t01 / cellCount << endl;
    cout << "    Seconds per non-zero expression matrix entry " << t01/double(nonZeroExpressionCount) << endl;
    cout << "    Seconds per inner loop iteration " << t01 / (double(nonZeroExpressionCount) * double(lshCount)) << endl;
    cout << "    Gflop/s " << 2. * 1e-9 * double(nonZeroExpressionCount) * double(lshCount) / t01 << endl;

}



// Compute the similarity (cosine of the angle) corresponding to each number of mismatching bits.
void Lsh::computeSimilarityTable()
{
    // Initialize the similarity table.
    const size_t lshCount = lshVectors.front().size();
    similarityTable.resize(lshCount + 1);

    // Loop over all possible numbers of mismatching bits.
    for(size_t mismatchingBitCount = 0;
        mismatchingBitCount <= lshCount; mismatchingBitCount++) {

        // Compute the angle between the vectors corresponding to
        // this number of mismatching bits.
        const double angle = double(mismatchingBitCount) *
            boost::math::double_constants::pi / double(lshCount);

        // The cosine of the angle is the similarity for
        // this number of mismatcning bits.
        CZI_ASSERT(mismatchingBitCount < similarityTable.size());
        similarityTable[mismatchingBitCount] = std::cos(angle);
    }

}


// Compute the LSH similarity between two cells,
// specified by their ids local to the cell set used by this Lsh object.
double Lsh::computeCellSimilarity(CellId localCellId0, CellId localCellId1)
{
    // Access the LSH signatures for the two cells.
    const BitSetInMemory signature0 = getSignature(localCellId0);
    const BitSetInMemory signature1 = getSignature(localCellId1);

    // Count the number of bits where the signatures of these two cells disagree.
    const size_t mismatchingBitCount = countMismatches(signatureWordCount, signature0, signature1);

    // Return the similarity corresponding to this number of mismatching bits.
    return similarityTable[mismatchingBitCount];
}



// Write to a csv file statistics of the cell LSH signatures..
void Lsh::writeSignatureStatistics(const string& csvFileName)
{
    ofstream csv(csvFileName);
    writeSignatureStatistics(csv);
}
void Lsh::writeSignatureStatistics(ostream& csv)
{
    csv << "Bit,Set,Unset,Total\n";

    for(size_t i = 0; i < info->lshCount; i++) {

        // Count the number of cells that have this bit set.
        size_t setCount = 0;
        for(CellId cellId = 0; cellId < info->cellCount; cellId++) {
            if(getSignature(cellId).get(i)) {
                ++setCount;
            }
        }
        const size_t unsetCount = info->cellCount - setCount;

        csv << i << "," << setCount << "," << unsetCount << "," << info->cellCount << "\n";

    }

}


void Lsh::remove()
{
    signatures.remove();
    info.remove();
}

