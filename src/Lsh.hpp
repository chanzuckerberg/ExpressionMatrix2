// Class to describe an RNA expression matrix.

#ifndef CZI_EXPRESSION_MATRIX2_LSH_HPP
#define CZI_EXPRESSION_MATRIX2_LSH_HPP


// Class Lsh is used to do quickly find pairs of similar cells
// using the cosine distance of Locality Sensitive Hashing (LSH)

#include "BitSet.hpp"
#include "Ids.hpp"
#include "MemoryMappedObject.hpp"
#include "MemoryMappedVector.hpp"

#include "cstddef.hpp"
#include "cstdint.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class ExpressionMatrixSubset;
        class Lsh;
    }
}



class ChanZuckerberg::ExpressionMatrix2::Lsh {
public:

    // Create a new Lsh object and store it on disk.
    // This can be expensive as it requires creating LSH signatures
    // for all cells in the specified cell set.
    Lsh(
        const string& name,             // Name prefix for memory mapped files.
        const ExpressionMatrixSubset&,  // For a subset of genes and cells.
        size_t lshCount,                // Number of LSH hyperplanes
        uint32_t seed                   // Seed to generate LSH hyperplanes.
        );

    // Access an existing Lsh object.
    Lsh(
        const string& name              // Name prefix for memory mapped files.
        );

    // Remove the memory mapped files.
    void remove();

    // Compute the LSH similarity between two cells,
    // specified by their ids local to the cell set used by this Lsh object.
    double computeCellSimilarity(CellId localCellId0, CellId localCellId1);
    size_t computeMismatchCount(CellId localCellId0, CellId localCellId1);


    // Get the signature corresponding to a given CellId (local to the cell set).
    BitSetPointer getSignature(CellId cellId)
    {
        const size_t offset = cellId*signatureWordCount;    // Offset of the signature of this cell (in 64 bit words)
        size_t* pointer = &(signatures[offset]);            // Pointer to the signature of this cell (as uint64_t words)
        return BitSetPointer(pointer, signatureWordCount);
    }

    void writeSignatureStatistics(const string& csvFileName);
    void writeSignatureStatistics(ostream&);

    CellId cellCount() const
    {
        return CellId(info->cellCount);
    }
    size_t lshCount() const
    {
        return info->lshCount;
    }
    size_t wordCount() const
    {
        return signatureWordCount;
    }

private:

    // The LSH vectors.
    // These are unit vectors in the Euclidean space of dimension equal to
    // the number of genes. Each defines a hyperplane orthogonal to it.
    // Indexed by [localGeneId][lshVectorId], where the localGeneId
    // is the index of the gene in the gene subset we are using,
    // and lshVectorId identifies the LSH hyperplanes and goes from 0 to lshCount.
    // This way, the components of the LSH vectors for a given gene are
    // all located contiguously in memory. This important to optimize
    // the speed of the computation of cell LSH signatures.
    vector< vector<double> > lshVectors;
    void generateLshVectors(
        size_t geneCount,
        size_t lshCount,                // Number of LSH hyperplanes
        uint32_t seed                   // Seed to generate LSH hyperplanes.
    );

    // The number of 64 bit words in each cell signature.
    size_t signatureWordCount;

    // The LSH signatures of all cells in the cell set we are using.
    // Each cell signature is a vector of bits.
    // The bit is 1 if the scalar product of the cell shifted expression vector
    // with the LSH vector corresponding to the bit position is positive,
    // and negative otherwise.
    MemoryMapped::Vector<uint64_t> signatures;
    void computeCellLshSignatures(const string& name, const ExpressionMatrixSubset&);

    // The similarity (cosine of the angle) corresponding to each number of mismatching bits.
    vector<double> similarityTable;
    void computeSimilarityTable();

    // Other information for this Lsh object.
    class Info {
    public:
        size_t cellCount;
        size_t lshCount;
    };
    MemoryMapped::Object<Info> info;

};

#endif
