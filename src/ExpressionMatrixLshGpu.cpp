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
    unsigned int seed               // The seed used to generate the LSH vectors.
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
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, k, geneSet, cellSet);

    // Vectors reused for each cell in the loop over cells below.
    vector<uint16_t> mismatchCounts(cellCount);
    vector< pair<uint16_t, CellId> > neighbors; // pair(mismatchCount, cellId1)

    // Compute the threshold on the number of mismatches.
    const size_t mismatchCountThreshold =
        lsh.computeMismatchCountThresholdFromSimilarityThreshold(similarityThreshold);



    // Loop over all cells.
    lsh.setupGpuKernel0();
    const auto t3 = std::chrono::steady_clock::now();
    for(CellId cellId0=0; cellId0!=cellCount; cellId0++) {

        // Use gpuKernel0 to compute the number of mismatches with all cells.
        lsh.gpuKernel0(cellId0, mismatchCounts);

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
    const auto t4 = std::chrono::steady_clock::now();
    const double t34 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count());
    cout << "Similar pairs computation on GPU took " << t34 << " s at ";
    cout << 1.e9*t34/(double(cellCount)*double(cellCount)) << " ns/pair" << endl;
    lsh.cleanupGpuKernel0();


    const auto t5 = std::chrono::steady_clock::now();
    const double t05 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t0)).count());
    cout << timestamp << "ExpressionMatrix::findSimilarPairs4Gpu ends. Took " << t05 << " s." << endl;
}

#endif
