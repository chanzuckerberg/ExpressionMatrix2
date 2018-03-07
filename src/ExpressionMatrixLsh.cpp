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
#include "charikar.hpp"
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

#include "algorithm.hpp"
#include <cmath>
#include "fstream.hpp"
#include <chrono>
#include <numeric>
#include <queue>
#include <random>



// Analyze the quality of a set of similar pairs.
void ExpressionMatrix::analyzeSimilarPairs(
    const string& similarPairsName,
    double csvDownsample) const
{
    // Open the SimilarPairs object we want to analyze.
    const SimilarPairs similarPairs(directoryName, similarPairsName, true);
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



// Used to test improvements to findSimilarPairs3.
void ExpressionMatrix::findSimilarPairs4(
    ostream& out,
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed               // The seed used to generate the LSH vectors.
    )
{
    out << timestamp << "ExpressionMatrix::findSimilarPairs4 begins." << endl;

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
    out << timestamp << "Creating expression matrix subset." << endl;
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
    out << timestamp << "Begin computing similarities for all cell pairs." << endl;
    const auto t0 = std::chrono::steady_clock::now();
    const CellId blockSize = 64;
    size_t pairCount = 0;
    size_t totalPairCount = size_t(cellCount)*(size_t(cellCount-1))/2;
    size_t blockCount = 0;
    for(CellId begin0=0; begin0<cellCount; begin0+=blockSize) {
        const CellId end0 = min(begin0+blockSize, cellCount);
        for(CellId begin1=0; begin1<=begin0; begin1+=blockSize) {
            if(blockCount>0 && ((blockCount%1000000)==0)) {
                out << timestamp << "Pair computation ";
                out << 100.*double(pairCount)/double(totalPairCount);
                out << "% complete." << endl;
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
    out << "Time for all pairs: " << t01 << " s." << endl;
    out << "Time per pair: " << t01/(0.5*double(cellCount)*double(cellCount-1)) << " s." << endl;

    // Store the pairs in a SimilarPairs object.
    out << timestamp << "Initializing SimilarPairs object." << endl;
    SimilarPairs similarPairs(directoryName, similarPairsName, geneSetName, cellSetName, k);
    out << timestamp << "Copying similar pairs." << endl;
    similarPairs.copy(tmp);


    // Sort the similar pairs for each cell by decreasing similarity.
    out << timestamp << "Sorting similar pairs." << endl;
    similarPairs.sort();
    out << timestamp << "ExpressionMatrix::findSimilarPairs4 ends." << endl;

    lsh.remove();

}
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
    findSimilarPairs4(cout, geneSetName, cellSetName, similarPairsName,
        k, similarityThreshold, lshCount, seed);
}


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
    SimilarPairs similarPairs(directoryName, similarPairsName, geneSetName, cellSetName, k);
    cout << timestamp << "Copying similar pairs." << endl;
    similarPairs.copy(tmp);


    // Sort the similar pairs for each cell by decreasing similarity.
    cout << timestamp << "Sorting similar pairs." << endl;
    similarPairs.sort();
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << timestamp << "ExpressionMatrix::findSimilarPairs5 ends. Took " << t01 << " s." << endl;

}



// Find similar cell pairs using LSH, without looping over all pairs.
// Like findSimilarPairs5, but using variable lsh slice length.
void ExpressionMatrix::findSimilarPairs7(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& lshName,          // The name of the Lsh object to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    const vector<int>& lshSliceLengths, // The number of bits in each LSH signature slice, in decreasing order.
    CellId maxCheck,                // Maximum number of cells to consider for each cell.
    size_t log2BucketCount
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs7 begins." << endl;
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
    if(lsh.cellCount() != cellCount) {
        throw runtime_error("LSH object " + lshName + " has a number of cells inconsistent with cell set " + cellSetName);
    }
    const size_t lshBitCount = lsh.lshCount();
    cout << "Number of LSH signature bits is " << lshBitCount << endl;

    // Check that the slice lengths are in decreasing order.
    const size_t sliceLengthCount = lshSliceLengths.size();
    for(size_t i=1; i<sliceLengthCount; i++) {
        if(lshSliceLengths[i] >= lshSliceLengths[i-1]) {
            throw runtime_error("The slice lengths are not in decreasing order.");
        }
    }

    // Check that the slice lengths are no more than 64.
    for(size_t i=0; i<sliceLengthCount; i++) {
        if(lshSliceLengths[i] > 64) {
            throw runtime_error("Each slice length can be at most 64 bits.");
        }
    }

    // Create SimilarPairs object that will store the results.
    SimilarPairs similarPairs(directoryName, similarPairsName, geneSetName, cellSetName, k);



    // For each slice length and signature slice of that length,
    // each cell is assigned to a bucket
    // based on the value of its signature slice.
    // The cells in each bucket will be stored in
    // tables4[sliceLengthId][sliceId][bucketId],
    // where:
    // - sliceLengthId corresponds to the slice lengths to be used,
    //   stored in lshSliceLengths in decreasing order.
    // - sliceId identifies the particular slice of that length
    //   (we have a total lshBitCount bits, which can be used
    //   to form lshBitCount/lshSliceLength possible signature
    //   slices each lshSliceLength bits in length).
    // - bucketId identifies the bucket.
    vector< vector< vector< vector<CellId> > > > table4;

    // Vector to contain the signature bits of each signature slice.
    // Indexed by [sliceLengthId][sliceId].
    vector< vector< vector< size_t > > > sliceBits3;

    // Assign cells to buckets.
    findSimilarPairs7AssignCellsToBuckets(lsh, lshSliceLengths, sliceBits3, table4, log2BucketCount);



    // Bit set to keep track which cellId1 cells we have already
    // looked at, for a given cellId0.
    BitSet cellMap(cellCount);

    // Other vectors used over and over again for each cell.
    vector<CellId> candidateNeighbors;
    vector< pair<uint32_t, CellId> > neighbors; // pair(mismatchCount, cellId1)

    const size_t mismatchCountThreshold =
        lsh.computeMismatchCountThresholdFromSimilarityThreshold(similarityThreshold);
    cout << "Mismatch count threshold is " << mismatchCountThreshold << endl;



    // For each cell, look at cells in the same bucket.
    // Stop when we found enough similar cells.
    cout << timestamp << "Finding similar cell pairs." << endl;
    const uint64_t bucketCount = (1ULL << log2BucketCount);
    const uint64_t bucketMask = bucketCount - 1ULL;
    for(CellId cellId0=0; cellId0<cellCount; cellId0++) {
        if(cellId0!=0 && (cellId0 % 1000)==0) {
            cout << timestamp << "Working on cell " << cellId0 << " of " << cellCount << endl;
        }
        const BitSetPointer signature = lsh.getSignature(cellId0);

        // Loop over slice lengths.
        for(size_t sliceLengthId=0; sliceLengthId<sliceLengthCount; sliceLengthId++) {
            const auto& table3 = table4[sliceLengthId];
            const auto& sliceBits2 = sliceBits3[sliceLengthId];

            // Extract the slice length.
            const size_t sliceLength = lshSliceLengths[sliceLengthId];

            // Compute the number of possible slices for this length.
            const size_t sliceCount = lshBitCount / sliceLength;

            // Loop over all possible signature slices of this length.
            for(size_t sliceId=0; sliceId<sliceCount; sliceId++) {
                const auto& table2 = table3[sliceId];
                const auto& sliceBits1 = sliceBits2[sliceId];

                // Extract this signature slice for this cell.
                const uint64_t signatureSlice = signature.getBits(sliceBits1);

                // Find the bucket that corresponds to this signature slice.
                const uint64_t bucketId =
                    (sliceLength<log2BucketCount) ?
                    signatureSlice :
                    (MurmurHash64A(&signatureSlice, 8, 231) & bucketMask);
                CZI_ASSERT(bucketId < table2.size());
                const auto& table1 = table2[bucketId];

                // Loop over cells in the same bucket.
                for(const CellId cellId1: table1) {
                    if(cellId1 == cellId0){
                        continue;
                    }
                    if(cellMap.get(cellId1)) {
                        continue;   // We already looked at this one.
                    }
                    cellMap.set(cellId1);
                    candidateNeighbors.push_back(cellId1);
                    const uint32_t mismatchCount = uint32_t(lsh.computeMismatchCount(cellId0, cellId1));
                    if(mismatchCount < mismatchCountThreshold) {
                        neighbors.push_back(make_pair(mismatchCount, cellId1));
                    }
                    if(candidateNeighbors.size() == maxCheck) {
                        break;
                    }
                }
                if(candidateNeighbors.size() == maxCheck) {
                    break;
                }
            }
            if(candidateNeighbors.size() == maxCheck) {
                break;
            }
        }

        // Only keep the k best neighbors, then sort them.
        // This is faster than sorting, then keeping the k best,
        // because it avoids doing a complete sorting of all of the neighbors.
        // Instead, keepBest uses std::nth_element, which does a partial sorting.
        keepBest(neighbors, k, std::less< pair<uint32_t, CellId> >());
        sort(neighbors.begin(), neighbors.end());

        // Store.
        for(const auto& neighbor: neighbors) {
            const CellId cellId1 = neighbor.second;
            const uint32_t mismatchCount = neighbor.first;
            const double similarity = lsh.getSimilarity(mismatchCount);
            similarPairs.addUnsymmetricNoCheck(cellId0, cellId1, similarity);
        }

        // Clean up our data structures so we can reuse them for the next cell.
        for(const CellId cellId1: candidateNeighbors) {
            cellMap.clear(cellId1);
        }
        candidateNeighbors.clear();
        neighbors.clear();

    }


    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << timestamp << "ExpressionMatrix::findSimilarPairs7 ends. Took " << t01 << " s." << endl;
}



// In the initial phase of findSimilarPairs7, we assign cells to buckets.
// For each slice length and signature slice of that length,
// each cell is assigned to a bucket
// based on the value of its signature slice.

// The cells in each bucket are stored in
// tables4[sliceLengthId][sliceId][bucketId],
// where:
// - sliceLengthId corresponds to the slice lengths to be used,
//   stored in lshSliceLengths in decreasing order.
// - sliceId identifies the particular slice of that length
//   (we have a total lshBitCount bits, which can be used
//   to form lshBitCount/lshSliceLength possible signature
//   slices each lshSliceLength bits in length).
// - bucketId identifies the bucket.

// Vector sliceBits3 is used to store the signature bits of each signature slice.
// Indexed by [sliceLengthId][sliceId].

void ExpressionMatrix::findSimilarPairs7AssignCellsToBuckets(
    Lsh& lsh,
    const vector<int>& lshSliceLengths,                     // The number of signature slice bits, in decreasing order.
    vector< vector< vector< size_t > > >& sliceBits3,       // The bits of each slice. See comments above.
    vector< vector< vector< vector<CellId> > > >& table4,   // The cells in each bucket. See comments above.
    size_t log2BucketCount
    )
{
    // Check that the slice lengths are in decreasing order.
    const size_t sliceLengthCount = lshSliceLengths.size();
    for(size_t i=1; i<sliceLengthCount; i++) {
        if(lshSliceLengths[i] >= lshSliceLengths[i-1]) {
            throw runtime_error("The slice lengths are not in decreasing order.");
        }
    }

    // Check that the slice lengths are no more than 64.
    for(size_t i=0; i<sliceLengthCount; i++) {
        if(lshSliceLengths[i] > 64) {
            throw runtime_error("Each slice length can be at most 64 bits.");
        }
    }


    // Initialize the table4 and sliceBits3 data structures.
    cout << timestamp << "Initializing data structures for findSimilarPairs7." << endl;
    table4.resize(sliceLengthCount);
    sliceBits3.resize(sliceLengthCount);
    const uint64_t bucketCount = (1ULL << log2BucketCount);
    const uint64_t bucketMask = bucketCount - 1ULL;
    const size_t lshBitCount = lsh.lshCount();
    for(size_t sliceLengthId=0; sliceLengthId<sliceLengthCount; sliceLengthId++) {
        auto& table3 = table4[sliceLengthId];
        auto& sliceBits2 = sliceBits3[sliceLengthId];

        // Extract the slice length.
        const size_t sliceLength = lshSliceLengths[sliceLengthId];

        // Compute the number of possible slices for this length.
        const size_t sliceCount = lshBitCount / sliceLength;
        table3.resize(sliceCount);
        sliceBits2.resize(sliceCount);
        const uint64_t tableSize = std::min(uint64_t(1ULL<<sliceLength), bucketCount);
        cout << "Number of slices of length " << sliceLength << " is " << sliceCount;
        cout << ". Table size is " << tableSize << endl;

        // Loop over all possible signature slices of this length.
        for(size_t sliceId=0; sliceId<sliceCount; sliceId++) {
            auto& table2 = table3[sliceId];
            table2.resize(tableSize);

            // Gather the signature bits.
            auto& sliceBits1 = sliceBits2[sliceId];
            sliceBits1.resize(sliceLength);
            size_t bitPosition = sliceId * sliceLength;
            for(size_t bitId=0; bitId<sliceLength; bitId++, ++bitPosition) {
                sliceBits1[bitId] = bitPosition;
            }
        }
    }



    // Assign cells to buckets.
    cout << timestamp << "Assigning cells to buckets." << endl;
    const CellId cellCount = lsh.cellCount();
    for(CellId cellId=0; cellId<cellCount; cellId++) {
        if(cellId!=0 && (cellId % 100000)==0) {
            cout << timestamp << "Working on cell " << cellId << " of " << cellCount << endl;
        }
        const BitSetPointer signature = lsh.getSignature(cellId);

        // Loop over slice lengths.
        for(size_t sliceLengthId=0; sliceLengthId<sliceLengthCount; sliceLengthId++) {
            auto& table3 = table4[sliceLengthId];
            const auto& sliceBits2 = sliceBits3[sliceLengthId];

            // Extract the slice length.
            const size_t sliceLength = lshSliceLengths[sliceLengthId];

            // Compute the number of possible slices for this length.
            const size_t sliceCount = lshBitCount / sliceLength;

            // Loop over all possible signature slices of this length.
            for(size_t sliceId=0; sliceId<sliceCount; sliceId++) {
                auto& table2 = table3[sliceId];
                const auto& sliceBits1 = sliceBits2[sliceId];

                // Extract this signature slice for this cell.
                const uint64_t signatureSlice = signature.getBits(sliceBits1);

                // Add this cell to the bucket that corresponds to this signature slice.
                const uint64_t bucketId =
                    (sliceLength<log2BucketCount) ?
                    signatureSlice :
                    (MurmurHash64A(&signatureSlice, 8, 231) & bucketMask);
                CZI_ASSERT(bucketId < table2.size());
                table2[bucketId].push_back(cellId);
            }
        }

    }
}



// Find similar cell pairs using LSH and the Charikar algorithm.
// See M. Charikar, "Similarity Estimation Techniques from Rounding Algorithms", 2002,
// section "5. Approximate Nearest neighbor Search in Hamming Space.".
// The Charikar algorithm is for approximate nearest neighbor, but with appropriate
// choices of the algorithm parameters permutationCount and searchCount
// can be used for approximate k nearest neighbors.
// In the Charikar paper, permutationCount is N and searchCount is 2N.
// To reduce memory requirements, we don't store all bits all bits of
// the permuted signatures - only the most significant permutedBitCount.
// In practice it is best to set this to 64, so the permuted signatured
// use only one 64-bit word each.
void ExpressionMatrix::findSimilarPairs6(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& lshName,          // The name of the Lsh object to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    size_t permutationCount,        // The number of bit permutations for the Charikar algorithm.
    size_t searchCount,             // The number of cells checked for each cell, in the Charikar algorithm.
    size_t permutedBitCount,        // The number of most significant bits stored for each permuted signature.
    int seed                        // The seed used to randomly generate the bit permutations.
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs6 begins." << endl;
    bool debug = false;
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

    // Sanity check on the number of most significant bits stored for
    // each permuted signature.
    if(permutedBitCount > lshCount) {
        throw runtime_error(
            "Argument permutationStoreBitCount " +
            to_string(permutedBitCount) +
            " exceeds number of signature bits " +
            to_string(lshCount));
    }


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
    // - The permuted signatures, in sorted order (first permutationStoreBitCount bits only).
    // - The corresponding cell ids, in order consistent with the permuted signatures.
    const size_t permutedWordCount = ((permutedBitCount-1) >> 6) + 1;
    cout << "Allocating " << ((8*permutationCount*size_t(lsh.cellCount())*permutedWordCount) >> 30) << " GB for permutation data." << endl;
    vector<Charikar::PermutationData> permutationData(
        permutationCount,
        Charikar::PermutationData(lsh.cellCount(), permutedWordCount));



    // For each of the permutations, compute permuted/sorted signatures.
    // We only compute and store the first permutedBitCount bits of each permuted signature.
    cout << timestamp << "Phase 1 of Charikar algorithm begins." << endl;
    const auto t1 = std::chrono::steady_clock::now();
    for(size_t permutationId=0; permutationId<permutationCount; permutationId++) {
        cout << timestamp << "Working on permutation " << permutationId << " of " << permutationCount << endl;

        // Generate a random permutation of the signature bits.
        vector<uint64_t> bitPermutation(lshCount);
        std::iota(bitPermutation.begin(), bitPermutation.end(), 0ULL);
        std::shuffle(bitPermutation.begin(), bitPermutation.end(), randomGenerator);
        bitPermutation.resize(permutedBitCount);    // Only keep the permutedBitCount most significant bits.

        if(debug) {
            cout << "Creating permutation data for permutation " << permutationId << ".\n";
            cout << "Bit permutation:\n";
            for(size_t i=0; i<lshCount; i++) {
                cout << i << " " << bitPermutation[i] << "\n";
            }
        }

        // Compute the permuted signatures for this permutation.
        BitSets permutedSignatures(cellCount, permutedWordCount);
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
        BitSets& thisPermutationBitSets = permutationData[permutationId].signatures;
        vector<CellId>& thisPermutationCellIds = permutationData[permutationId].cellIds;
        CZI_ASSERT(thisPermutationBitSets.bitSetCount == cellCount);
        CZI_ASSERT(thisPermutationBitSets.wordCount == permutedWordCount);
        CZI_ASSERT(thisPermutationCellIds.size() == cellCount);
        for(CellId i=0; i<cellCount; i++) {
            pair<BitSetPointer, CellId>& p = table[i];
            thisPermutationBitSets.set(i, p.first);
            thisPermutationCellIds[i] = p.second;
        }
        permutationData[permutationId].computeCellPositions();

        if(debug) {
            cout << "Permutation data for permutation " << permutationId << ":" << endl;
            for(CellId i=0; i<cellCount; i++) {
                cout << thisPermutationBitSets[i].getString(lshCount) << " " << thisPermutationCellIds[i] << "\n";
            }
        }
    }
    const auto t2 = std::chrono::steady_clock::now();
    cout << timestamp << "Phase 1 of Charikar algorithm took ";
    cout << 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count()) << " s." << endl;

    // Temporary storage of pairs for each cell.
    vector< vector< pair<CellId, float> > > pairs(cellCount);
    vector< pair<CellId, float> > cellNeighbors;


    // At this point the necessary data structures are in place and we can use the Charikar algorithm
    // to find the neighbors of each cell.
    debug = false;
    cout << timestamp << "Phase 2 of Charikar algorithm begins." << endl;
    const auto t3 = std::chrono::steady_clock::now();
    for(CellId cellId0=0; cellId0<cellCount; cellId0++) {
        if(cellId0!=0 && (cellId0%100000)==0) {
            cout << timestamp << "Working on cell " << cellId0 << " of " << cellCount << endl;
        }

        // Extract the permuted signatures for this cell.
        vector<BitSetPointer> signatures0(permutationCount);
        for(size_t permutationId=0; permutationId<permutationCount; permutationId++) {
           Charikar::PermutationData& p = permutationData[permutationId];
           const size_t i = p.cellPositions[cellId0];
           CZI_ASSERT(p.cellIds[i] == cellId0);
           signatures0[permutationId] = p.signatures[i];
        }


        // Create the priority queue of Charikar pointers.
        std::priority_queue<Charikar::Pointer> priorityQueue;
        for(size_t permutationId=0; permutationId<permutationCount; permutationId++) {
            Charikar::PermutationData& pd = permutationData[permutationId];
            const size_t i = pd.cellPositions[cellId0];
            CZI_ASSERT(pd.cellIds[i] == cellId0);

            // Add the forward moving pointer.
            if(i < cellCount-1) {
                Charikar::Pointer pointer;
                pointer.permutationId = permutationId;
                pointer.index = i+1;
                pointer.movesForward = true;
                pointer.prefixLength = commonPrefixLength(
                    signatures0[permutationId], pd.signatures[i+1]);
                priorityQueue.push(pointer);
            }

            // Add the backward moving pointer.
            if(i > 1) {
                Charikar::Pointer pointer;
                pointer.permutationId = permutationId;
                pointer.index = i-1;
                pointer.movesForward = false;
                pointer.prefixLength = commonPrefixLength(
                    signatures0[permutationId], pd.signatures[i-1]);
                priorityQueue.push(pointer);
            }

        }



        // The heart of the Charikar algorithm begins here.
        // At each iteration we get the pointer with the best prefix.
        cellNeighbors.clear();
        for(size_t iteration=0; iteration<searchCount; iteration++) {

            // Get the pointer with the best prefix.
            if(priorityQueue.empty()) {
                break;
            }
            Charikar::Pointer pointer = priorityQueue.top();
            priorityQueue.pop();

            // Compute the number of mismatches.
            Charikar::PermutationData& pd = permutationData[pointer.permutationId];
            const CellId cellId1 = pd.cellIds[pointer.index];
            CZI_ASSERT(cellId1 != cellId0); // By construction.
            const double similarity = lsh.computeCellSimilarity(cellId0, cellId1);
            if(similarity > similarityThreshold) {
                cellNeighbors.push_back(make_pair(cellId1, similarity));
                if(false) {
                    cout << cellId0 << " " << cellId1 << " " << pointer.prefixLength << " " << similarity << endl;
                }
            }

            // Update the pointer and requeue it.
            if(pointer.movesForward) {
                if(pointer.index < cellCount-1) {
                    ++pointer.index;
                    pointer.prefixLength = commonPrefixLength(
                        signatures0[pointer.permutationId], pd.signatures[pointer.index]);
                    priorityQueue.push(pointer);
                }
            } else {
                if(pointer.index > 0) {
                    --pointer.index;
                    pointer.prefixLength = commonPrefixLength(
                        signatures0[pointer.permutationId], pd.signatures[pointer.index]);
                    priorityQueue.push(pointer);
                }
            }

        }

        // Sort, deduplicate, keep the best k.
        sort(cellNeighbors.begin(), cellNeighbors.end(),
            OrderPairsBySecondGreaterThenByFirstLess< pair<CellId, float> >());
        cellNeighbors.resize(unique(cellNeighbors.begin(), cellNeighbors.end()) - cellNeighbors.begin());
        if(cellNeighbors.size() > k) {
            cellNeighbors.resize(k);
        }

        // Store.
        pairs[cellId0] = cellNeighbors;
    }
    const auto t4 = std::chrono::steady_clock::now();
    cout << timestamp << "Phase 2 of Charikar algorithm took ";
    cout << 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count()) << " s." << endl;



    // Store the pairs in a SimilarPairs object.
    cout << timestamp << "Initializing SimilarPairs object." << endl;
    SimilarPairs similarPairs(directoryName, similarPairsName, geneSetName, cellSetName, k);
    cout << timestamp << "Copying similar pairs." << endl;
    similarPairs.copy(pairs);


    // Sort the similar pairs for each cell by decreasing similarity.
    // (This should not be necessary - consider removing).
    cout << timestamp << "Sorting similar pairs." << endl;
    similarPairs.sort();

    const auto t5 = std::chrono::steady_clock::now();
    const double t05 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t0)).count());
    cout << timestamp << "ExpressionMatrix::findSimilarPairs6 ends. Took " << t05 << " s." << endl;
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
    const SimilarPairs similarPairs0(directoryName, similarPairsName0, true);
    const SimilarPairs similarPairs1(directoryName, similarPairsName1, true);

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

    // Gather cells with the same signature.
#if 0
    vector< pair<BitSetPointer, CellId> > table;
    table.reserve(cellCount);
    for(CellId cellId=0; cellId<cellCount; cellId++) {
        table.push_back(make_pair(lsh.getSignature(cellId), cellId));
    }
    cout << timestamp << "Sorting by signature." << endl;
    sort(table.begin(), table.end());
#endif
    map<BitSetPointer, vector<CellId> > signatureMap;
    for(CellId cellId=0; cellId<cellCount; cellId++) {
        signatureMap[lsh.getSignature(cellId)].push_back(cellId);
    }


    // Create a table of signatures ordered by decreasing number of cells.
    vector< pair<BitSetPointer, size_t> > signatureTable;
    for(const auto& p: signatureMap) {
        const BitSetPointer signature = p.first;
        const size_t size = p.second.size();
        signatureTable.push_back(make_pair(signature, size));
    }
    sort(signatureTable.begin(), signatureTable.end(),
        OrderPairsBySecondGreater< pair<BitSet, size_t> >());



    // Write a csv file with one line for each distinct signature.
    {
        ofstream csvOut("Signatures.csv");
        for(const auto& p: signatureTable) {
            const BitSetPointer signature = p.first;
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
