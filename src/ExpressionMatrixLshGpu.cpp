#if CZI_EXPRESSION_MATRIX2_BUILD_FOR_GPU


#include "ExpressionMatrix.hpp"
#include "ExpressionMatrixSubset.hpp"
#include "heap.hpp"
#include "Lsh.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



void ExpressionMatrix::findSimilarPairs4Gpu(
    const string& geneSetName,      // The name of the gene set to be used.
    const string& cellSetName,      // The name of the cell set to be used.
    const string& lshName,          // The name of the Lsh object to be used.
    const string& similarPairsName, // The name of the SimilarPairs object to be created.
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    size_t lshCount,                // The number of LSH vectors to use.
    unsigned int seed,              // The seed used to generate the LSH vectors.
    unsigned int kernel,            // The GPU kernel (algorithm) to use.
    CellId blockSize                // The number of cells processed by each kernel instance.
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs4Gpu begins." << endl;
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

    // Create the expression matrix subset for this gene set and cell set.
    cout << timestamp << "Creating expression matrix subset." << endl;
    const string expressionMatrixSubsetName =
        directoryName + "/tmp-ExpressionMatrixSubset-" + similarPairsName;
    ExpressionMatrixSubset expressionMatrixSubset(
        expressionMatrixSubsetName, geneSet, cellSet, cellExpressionCounts);

    // Create the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/Lsh-" + lshName);
    if(lsh.cellCount() != cellSet.size()) {
        throw runtime_error("LSH object " + lshName + " has a number of cells inconsistent with cell set " + cellSetName);
    }
    lsh.initializeGpu();
    cout << "GPU computing will use " << lsh.getGpuName() << "." << endl;

    // Load the LSH signatures to the GPU.
    const auto t1 = std::chrono::steady_clock::now();
    lsh.loadSignaturesToGpu();
    const auto t2 = std::chrono::steady_clock::now();
    const double t12 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count());
    cout << "Loading cell signatures to GPU took " << t12 << " s at ";
    cout << (double(cellCount)*double(lsh.wordCount())*sizeof(uint64_t) /
        (t12*1024.*1024.)) << " MB/s." << endl;

    // Create the SimilarPairs object that will store the results.
    SimilarPairs similarPairs(directoryName, similarPairsName, geneSetName, cellSetName, k);

    // Call the appropriate lower level function.
    switch(kernel) {
    case 0:
        findSimilarPairs4GpuKernel0(k, similarityThreshold, lsh, similarPairs);
        break;
    case 1:
        findSimilarPairs4GpuKernel1(k, similarityThreshold, lsh, similarPairs, blockSize);
        break;
    case 2:
        findSimilarPairs4GpuKernel2(k, similarityThreshold, lsh, similarPairs, blockSize);
        break;
    default:
        throw runtime_error(
            "Invalid kernel " + std::to_string(kernel) +
            " specified in findSimilarPairs4Gpu call. ");
    }

    const auto t3 = std::chrono::steady_clock::now();
    const double t03 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t0)).count());
    cout << timestamp << "ExpressionMatrix::findSimilarPairs4Gpu ends. Took " << t03 << " s." << endl;
}


// Kernel 0: compute the number of mismatches between the signature
// of cell cellId0 and the signatures of all other cells.
// This is a simple but working kernel.
// It has high overhead and low performance
// because of the small amount of work done by each instance.
void ExpressionMatrix::findSimilarPairs4GpuKernel0(
    size_t k,
    double similarityThreshold,
    Lsh& lsh,
    SimilarPairs& similarPairs)
{
    const CellId cellCount = lsh.cellCount();

    // Vectors reused for each cell in the loop over cells below.
    vector<uint16_t> mismatchCounts(cellCount);
    vector< pair<uint16_t, CellId> > neighbors; // pair(mismatchCount, cellId1)

    // Compute the threshold on the number of mismatches.
    const size_t mismatchCountThreshold =
        lsh.computeMismatchCountThresholdFromSimilarityThreshold(similarityThreshold);


    // Loop over all cells.
    lsh.setupGpuKernel0(mismatchCounts);
    const auto t0 = std::chrono::steady_clock::now();
    for(CellId cellId0=0; cellId0!=cellCount; cellId0++) {

        // Use gpuKernel0 to compute the number of mismatches with all cells.
        lsh.gpuKernel0(cellId0);

        // Candidate neighbors are the ones where the number of mismatches is small.
        neighbors.clear();
        for(CellId cellId1=0; cellId1!=cellCount; cellId1++) {
            if(cellId1 == cellId0) {
                continue;
            }
            const uint16_t mismatchCount = mismatchCounts[cellId1];
            if(mismatchCount <= mismatchCountThreshold) {
                neighbors.push_back(make_pair(mismatchCount, cellId1));
            }
        }

        // Only keep the k best neighbors, then sort them.
        // This is faster than sorting, then keeping the k best,
        // because it avoids doing a complete sorting of all of the neighbors.
        // Instead, keepBest uses std::nth_element, which does a partial sorting.
        keepBest(neighbors, k, std::less< pair<uint16_t, CellId> >());
        sort(neighbors.begin(), neighbors.end());

        // Store.
        for(const auto& neighbor: neighbors) {
            const CellId cellId1 = neighbor.second;
            const uint16_t mismatchCount = neighbor.first;
            const double similarity = lsh.getSimilarity(mismatchCount);
            similarPairs.addUnsymmetricNoCheck(cellId0, cellId1, similarity);
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << "Similar pairs computation on GPU took " << t01 << " s at ";
    cout << 1.e9*t01/(double(cellCount)*double(cellCount)) << " ns/pair" << endl;
    lsh.cleanupGpuKernel0();
}



// Kernel 1: compute the number of mismatches between the signature
// of cells in range [cellId0Begin, cellId0End]
// and the signatures of all other cells.
// This is parallelized over cellId1.
// This has lower overhead than kernel 0
// because each kernel instance does more work
// (processes n0 values of cellId0 instead of just one).
// For efficient memory access, the blockSize mismatch counts
// for each cellId1 are stored contiguously.
void ExpressionMatrix::findSimilarPairs4GpuKernel1(
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    Lsh& lsh,
    SimilarPairs& similarPairs,
    CellId blockSize)
{
    const CellId cellCount = lsh.cellCount();

    // Vectors reused for each cell in the loop over cells below.
    vector<uint16_t> mismatchCounts(cellCount * blockSize);
    vector< pair<uint16_t, CellId> > neighbors; // pair(mismatchCount, cellId1)

    // Compute the threshold on the number of mismatches.
    const size_t mismatchCountThreshold =
        lsh.computeMismatchCountThresholdFromSimilarityThreshold(similarityThreshold);


    // Loop over all cells.
    lsh.setupGpuKernel1(mismatchCounts, blockSize);
    const auto t0 = std::chrono::steady_clock::now();
    uint64_t blockId = 0;
    for(CellId cellId0Begin=0; cellId0Begin<cellCount; cellId0Begin+=blockSize, ++blockId) {
        if((blockId%1000) == 0) {
            cout << timestamp << " Working on cell " << cellId0Begin << " of " << cellCount << endl;
        }
        const CellId cellId0End = min(cellCount, cellId0Begin + blockSize);

        // Use gpuKernel0 to compute the number of mismatches
        // of cells in range [cellId0Begin, cellId0End) with all cells.
        // For efficient memory access in the GPU, the mismatch counts
        // for each cellId1 are stored contiguously,
        // with stride blockSize between consecutive values of cellId1.
        // const auto tA = std::chrono::steady_clock::now();
        lsh.gpuKernel1(cellId0Begin, cellId0End);
        // const auto tB = std::chrono::steady_clock::now();

        // Process the results for all the cellId0 values in this block.
        for(CellId cellId0=cellId0Begin; cellId0!=cellId0End; ++cellId0) {
            // const auto tt0 = std::chrono::steady_clock::now();

            // Candidate neighbors are the ones where the number of mismatches is small.
            neighbors.clear();
            // const auto tt1 = std::chrono::steady_clock::now();
            for(CellId cellId1=0; cellId1!=cellCount; cellId1++) {
                if(cellId1 == cellId0) {
                    continue;
                }
                const uint16_t mismatchCount = mismatchCounts[cellId1*blockSize + (cellId0-cellId0Begin)];
                if(mismatchCount <= mismatchCountThreshold) {
                    neighbors.push_back(make_pair(mismatchCount, cellId1));
                }
            }
            // const auto tt2 = std::chrono::steady_clock::now();

            // Only keep the k best neighbors, then sort them.
            // This is faster than sorting, then keeping the k best,
            // because it avoids doing a complete sorting of all of the neighbors.
            // Instead, keepBest uses std::nth_element, which does a partial sorting.
            keepBest(neighbors, k, std::less< pair<uint16_t, CellId> >());
            // const auto tt3 = std::chrono::steady_clock::now();
            sort(neighbors.begin(), neighbors.end());
            // const auto tt4 = std::chrono::steady_clock::now();

            // Store.
            for(const auto& neighbor: neighbors) {
                const CellId cellId1 = neighbor.second;
                const uint16_t mismatchCount = neighbor.first;
                const double similarity = lsh.getSimilarity(mismatchCount);
                similarPairs.addUnsymmetricNoCheck(cellId0, cellId1, similarity);
            }
            /*
            const auto tt5 = std::chrono::steady_clock::now();
            cout << (tt1-tt0).count() << " ";
            cout << (tt2-tt1).count() << " ";
            cout << (tt3-tt2).count() << " ";
            cout << (tt4-tt3).count() << " ";
            cout << (tt5-tt4).count() << "\n";
            */
        }
        // const auto tC = std::chrono::steady_clock::now();
        // cout << (tB-tA).count() << " " << (tC-tB).count() << "\n";
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << "Similar pairs computation on GPU took " << t01 << " s at ";
    cout << 1.e9*t01/(double(cellCount)*double(cellCount)) << " ns/pair" << endl;
    lsh.cleanupGpuKernel1();
}


void ExpressionMatrix::findSimilarPairs4GpuKernel2(
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    Lsh& lsh,
    SimilarPairs& similarPairs,
    CellId blockSize)
{
    const CellId cellCount = lsh.cellCount();

    // Compute the threshold on the number of mismatches.
    const size_t mismatchCountThreshold =
        lsh.computeMismatchCountThresholdFromSimilarityThreshold(similarityThreshold);

    // Vectors reused for each block in the loop below.
    vector<CellId> neighbors(k * (mismatchCountThreshold+1) * blockSize);
    vector<CellId> neighborCounts((mismatchCountThreshold+1) * blockSize);

    // Loop over blocks of cells.
    lsh.setupGpuKernel2(neighbors, neighborCounts, blockSize, k, mismatchCountThreshold);
    const auto t0 = std::chrono::steady_clock::now();
    uint64_t blockId = 0;
    for(CellId cellId0Begin=0; cellId0Begin<cellCount; cellId0Begin+=blockSize, ++blockId) {
        if(true /*(blockId%1000) == 0*/) {
            cout << timestamp << " Working on cell " << cellId0Begin << " of " << cellCount << endl;
        }
        const CellId cellId0End = min(cellCount, cellId0Begin + blockSize);

        // Run GPU kernel 2 for this block.
        const auto tA = std::chrono::steady_clock::now();
        lsh.gpuKernel2(cellId0Begin, cellId0End);
        const auto tB = std::chrono::steady_clock::now();



        // Process the results of this block.
        // Loop over cells in this block.
        for(CellId cellId0=cellId0Begin; cellId0!=cellId0End; cellId0++) {

            // Pointers to the portions of the neighbors and neighborCounts
            // buffer that belong to this cell.
            const uint32_t* neighbors0 = neighbors.data() +
                k * (mismatchCountThreshold+1) * (cellId0-cellId0Begin);
            const uint32_t* neighborCounts0 = neighborCounts.data() +
                (mismatchCountThreshold+1) * (cellId0-cellId0Begin);

            // TGotal number of neighbors we already stored for this cell,
            // for all mismatch counts.
            size_t totalNeighborCount = 0;

            // Loop over mismatch counts.
            for(size_t mismatchCount=0; mismatchCount<=mismatchCountThreshold;
                mismatchCount++, neighbors0+=k, neighborCounts0++) {
                const double similarity = lsh.getSimilarity(mismatchCount);

                // Loop over neighbors for this mismatch count.
                for(size_t i=0; i<*neighborCounts0; i++) {
                    const CellId cellId1 = neighbors0[i];
                    if(cellId1 == cellId0) {
                        continue;
                    }
                    similarPairs.addUnsymmetricNoCheck(cellId0, cellId1, similarity);
                    ++totalNeighborCount;
                    if(totalNeighborCount == k) {
                        break;
                    }
                }
                if(totalNeighborCount == k) {
                    break;
                }
            }
        }
        const auto tC = std::chrono::steady_clock::now();
        cout << (tB-tA).count() << " " << (tC-tB).count() << endl;
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << "Similar pairs computation on GPU took " << t01 << " s at ";
    cout << 1.e9*t01/(double(cellCount)*double(cellCount)) << " ns/pair" << endl;
    lsh.cleanupGpuKernel2();

}



// GPU version of findSimilarPairs7. See Lsh.cl for details.
void ExpressionMatrix::findSimilarPairs7Gpu(
    const string& geneSetName,          // The name of the gene set to be used.
    const string& cellSetName,          // The name of the cell set to be used.
    const string& lshName,              // The name of the Lsh object to be used.
    const string& similarPairsName,     // The name of the SimilarPairs object to be created.
    size_t k,                           // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,         // The minimum similarity for a pair to be stored.
    const vector<int>& lshSliceLengths, // The number of bits in each LSH signature slice, in decreasing order.
    CellId maxCheck,                    // Maximum number of candidate neighbors to consider for each cell.
    size_t log2BucketCount,
    CellId blockSize                    // The number of cells processed by each kernel instance on the GPU.
    )
{
    cout << timestamp << "ExpressionMatrix::findSimilarPairs7Gpu begins." << endl;
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

    // Create the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/Lsh-" + lshName);
    if(lsh.cellCount() != cellSet.size()) {
        throw runtime_error("LSH object " + lshName + " has a number of cells inconsistent with cell set " + cellSetName);
    }
    lsh.initializeGpu();
    cout << "GPU computing will use " << lsh.getGpuName() << "." << endl;

    // Load the LSH signatures to the GPU.
    const auto t1 = std::chrono::steady_clock::now();
    lsh.loadSignaturesToGpu();
    const auto t2 = std::chrono::steady_clock::now();
    const double t12 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count());
    cout << "Loading cell signatures to GPU took " << t12 << " s at ";
    cout << (double(cellCount)*double(lsh.wordCount())*sizeof(uint64_t) /
        (t12*1024.*1024.)) << " MB/s." << endl;

    // Create the SimilarPairs object that will store the results.
    SimilarPairs similarPairs(directoryName, similarPairsName, geneSetName, cellSetName, k);

    // Call the appropriate lower level function.
    // For now we only have kernel 3, but this allows for
    // alternate kernels, to facilitate experimenting.
    const int kernel = 3;
    switch(kernel) {
    case 3:
        findSimilarPairs7GpuKernel3(
            k, similarityThreshold,
            lsh, similarPairs, lshSliceLengths,
            maxCheck, size_t(blockSize), log2BucketCount);
        break;
    default:
        throw runtime_error(
            "Invalid kernel " + std::to_string(kernel) +
            " specified in findSimilarPairs7Gpu call. ");
    }

    const auto t3 = std::chrono::steady_clock::now();
    const double t03 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t0)).count());
    cout << timestamp << "ExpressionMatrix::findSimilarPairs7Gpu ends. Took " << t03 << " s." << endl;
}



void ExpressionMatrix::findSimilarPairs7GpuKernel3(
    size_t k,                       // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold,     // The minimum similarity for a pair to be stored.
    Lsh& lsh,
    SimilarPairs& similarPairs,
    const vector<int>& lshSliceLengths, // The number of bits in each LSH signature slice, in decreasing order.
    CellId maxCheck,
    CellId blockSize,
    size_t log2BucketCount)
{
    const CellId cellCount = lsh.cellCount();
    const size_t lshBitCount = lsh.lshCount();
    const size_t sliceLengthCount = lshSliceLengths.size();
    const uint64_t bucketCount = (1ULL << log2BucketCount);
    const uint64_t bucketMask = bucketCount - 1ULL;

    // Compute the threshold on the number of mismatches.
    const size_t mismatchCountThreshold =
        lsh.computeMismatchCountThresholdFromSimilarityThreshold(similarityThreshold);



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

    // Vectors reused for each block in the loop below.
    // Each of these corresponds to a buffer in the GPU.
    vector<CellId> neighbors(k * (mismatchCountThreshold+1) * blockSize);
    vector<CellId> neighborCounts((mismatchCountThreshold+1) * blockSize);
    vector<CellId> cellId1s(maxCheck * blockSize);
    vector<CellId> cellId1Counts(blockSize);

    // Set up corresponding buffers in the GPU.
    cout << timestamp << "Setting up buffers in the GPU." << endl;
    lsh.setupGpuKernel3(
        k, mismatchCountThreshold,
        maxCheck, blockSize,
        cellId1s, cellId1Counts,
        neighbors, neighborCounts);

    // Bit set to keep track which cellId1 cells we have already
    // added as a candidate neighbor, for a given cellId0.
    BitSet cellMap(cellCount);

    // Loop over blocks of cells.
    cout << timestamp << "Similar pairs computation begins." << endl;
    const auto t0 = std::chrono::steady_clock::now();
    uint64_t blockId = 0;
    for(CellId cellId0Begin=0; cellId0Begin<cellCount; cellId0Begin+=blockSize, ++blockId) {
        const auto tA = std::chrono::steady_clock::now();
        if(true /*(blockId%1000) == 0*/) {
            cout << timestamp << " Working on cell " << cellId0Begin << " of " << cellCount << endl;
        }
        const CellId cellId0End = min(cellCount, cellId0Begin + blockSize);



        // Find the candidate neighbor cellId1s for each value of cellId0 in this block.
        // For each cell, look at cells in the same bucket.
        for(CellId cellId0=cellId0Begin; cellId0!=cellId0End; cellId0++) {
            const size_t i0 = cellId0 - cellId0Begin;
            CellId count0 = 0;
            CellId* const candidates0 = cellId1s.data() + maxCheck * i0;

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
                    // const auto tt0 = std::chrono::steady_clock::now();
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
                            continue;   // We already added this one.
                        }
                        cellMap.set(cellId1);
                        candidates0[count0++] = cellId1;
                        if(count0 == maxCheck) {
                            break;
                        }
                    }
                    if(count0 == maxCheck) {
                        break;
                    }
                }
                if(count0 == maxCheck) {
                    break;
                }
            }


            // Store the number of candidate neighbors.
            cellId1Counts[i0] = count0;

            // Clean up our bit set so we can reuse it for the next cell.
            for(size_t i1=0; i1<count0; i1++) {
                cellMap.clear(candidates0[i1]);
            }
        }



        // Run GPU kernel 3 for this block.
        const auto tB = std::chrono::steady_clock::now();
        lsh.gpuKernel3(cellId0Begin, cellId0End);
        const auto tC = std::chrono::steady_clock::now();



        // Process the results of this block.
        // Loop over cells in this block.
        for(CellId cellId0=cellId0Begin; cellId0!=cellId0End; cellId0++) {

            // Pointers to the portions of the neighbors and neighborCounts
            // buffer that belong to this cell.
            const uint32_t* neighbors0 = neighbors.data() +
                k * (mismatchCountThreshold+1) * (cellId0-cellId0Begin);
            const uint32_t* neighborCounts0 = neighborCounts.data() +
                (mismatchCountThreshold+1) * (cellId0-cellId0Begin);

            // Total number of neighbors we already stored for this cell,
            // for all mismatch counts.
            size_t totalNeighborCount = 0;

            // Loop over mismatch counts.
            for(size_t mismatchCount=0; mismatchCount<=mismatchCountThreshold;
                mismatchCount++, neighbors0+=k, neighborCounts0++) {
                const double similarity = lsh.getSimilarity(mismatchCount);

                // Loop over neighbors for this mismatch count.
                for(size_t i=0; i<*neighborCounts0; i++) {
                    const CellId cellId1 = neighbors0[i];
                    if(cellId1 == cellId0) {
                        continue;
                    }
                    similarPairs.addUnsymmetricNoCheck(cellId0, cellId1, similarity);
                    ++totalNeighborCount;
                    if(totalNeighborCount == k) {
                        break;
                    }
                }
                if(totalNeighborCount == k) {
                    break;
                }
            }
        }
        const auto tD = std::chrono::steady_clock::now();
        cout << (tB-tA).count() << " " << (tC-tB).count() << " " << (tD-tC).count() << endl;
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
    cout << timestamp << "Similar pairs computation on GPU took " << t01 << " s at ";
    cout << 1.e9*t01/(double(cellCount)*double(cellCount)) << " ns/pair" << endl;
    lsh.cleanupGpuKernel3();
}
#endif
