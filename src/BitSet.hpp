#ifndef CZI_EXPRESSION_MATRIX2_BIT_SET_HPP
#define CZI_EXPRESSION_MATRIX2_BIT_SET_HPP


// A bare bone bit set class with just the functionality needed for operations
// on LSH signatures. Similar to boost::dynamic_bitset.

#include "CZI_ASSERT.hpp"

#include "algorithm.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class BitSetPointer;
        class BitSet;
        uint64_t countMismatches(const BitSetPointer&, const BitSetPointer&);
    }
}



// Class that stores pointers to the begin and end of a bit set.
// Does not own the memory, and copies are shallow copies.
// Implements low level operations on bit sets.
class ChanZuckerberg::ExpressionMatrix2::BitSetPointer {
public:

    // Begin and and pointers of the bit set.
    uint64_t* begin;
    uint64_t* end;

    // Constructors.
    BitSetPointer(uint64_t* begin=0, uint64_t* end=0) : begin(begin), end(end) {}
    BitSetPointer(uint64_t* begin, uint64_t wordCount) : begin(begin), end(begin+wordCount) {}

    // Return the number of 64-bit words in the bit set.
    uint64_t wordCount() const
    {
        return end - begin;
    }

    // Get the index of the word that contains the bit
    // corresponding to a given bit position.
    uint64_t getWordIndex(uint64_t bitPosition) const
    {
        const uint64_t wordIndex = bitPosition >> 6;
        CZI_ASSERT(wordIndex < wordCount()); // Skip this check, for speed.
        return wordIndex;
    }

    // Get the position in its word of the bit corresponding to a given bit position.
    uint64_t getBitPosition(uint64_t bitPosition) const
    {
        // The first bit is in the most significant position,
        // so sorting the bit set using integer comparison
        // does a lexicographical ordering of the bit set.
        return 63ULL - (bitPosition & 63ULL);
    }

    // Get the bit at a given position.
    bool get(uint64_t bitPosition) const
    {
        // Find the word containing the bit.
        const uint64_t wordIndex = getWordIndex(bitPosition);
        const uint64_t& word = begin[wordIndex];

        // Find the position of this bit in the word.
        const uint64_t bitPositionInWord = getBitPosition(bitPosition);

        // Test the bit.
        const uint64_t mask = 1ULL << bitPositionInWord;
        return (word & mask) != 0ULL;
    }

    // Set the bit at a given position.
    void set(uint64_t bitPosition)
    {
        // Find the word containing the bit.
        const uint64_t wordIndex = getWordIndex(bitPosition);
        uint64_t& word = begin[wordIndex];

        // Find the position of this bit in the word.
        const uint64_t bitPositionInWord = getBitPosition(bitPosition);

        // Set the bit.
        word |= 1ULL << bitPositionInWord;
    }

    // Get a uint64_t containing bits at specified positions.
    // The bits are specified in a vector.
    // The last specified bit goes in the least significant position
    // of the return value.
    uint64_t getBits(const vector<size_t>& bitPositions) const
    {
        uint64_t bits = 0;
        for(const size_t bitPosition: bitPositions) {
            bits <<= 1;
            bits += get(bitPosition);
        }
        return bits;
    }

    // Get a string representing a specified number of bits of the bit set.
    // The string contains an underscore for each zero bit
    // and an x for each one bit.
    string getString(uint64_t bitCount) const
    {
        string s;
        for(uint64_t i=0; i<bitCount; i++) {
            if(get(i)) {
                s += 'x';
            } else {
                s += '_';
            }
        }
        return s;

    }

    // Fill the bits of this bit set using a permutation
    // of the bits of another bit set.
    void fillUsingPermutation(
        const vector<int>& permutation,
        BitSetPointer& that)
    {
        clear();
        for(size_t i=0; i<permutation.size(); i++) {
            if(that.get(permutation[i])) {
                set(i);
            }
        }
    }

    void clear()
    {
        fill(begin, end, 0ULL);
    }

    bool operator<(const BitSetPointer& that) const
    {
        return std::lexicographical_compare(begin, end, that.begin, that.end);
    }
};



// Bit set that owns its memory.
class ChanZuckerberg::ExpressionMatrix2::BitSet : public BitSetPointer {
public:


    // Create a bit set with all bits set to zero.
    BitSet(uint64_t bitCount)
    {
        CZI_ASSERT(bitCount > 0);
        const uint64_t wordCount = ((bitCount - 1ULL) >> 6ULL) + 1ULL;
        begin = new uint64_t(wordCount);
        end = begin + wordCount;
        clear();
    }

    // Create a bit set as a copy of the bit set described by a BitSetPointer.
    BitSet(const BitSetPointer& that)
    {
        const uint64_t wordCount = that.wordCount();
        begin = new uint64_t(wordCount);
        end = begin + wordCount;
        copy(that.begin, that.end, begin);
    }

    // The BitSet owns its memory, so we must provide
    // a destructor, copy constructor, and assignment operator.
    ~BitSet()
    {
        delete[] begin;
    }
    BitSet(const BitSet& that) : BitSet(BitSetPointer(that)) {}
    BitSet& operator=(const BitSet& that)
    {
        delete[] begin;
        const uint64_t wordCount = that.wordCount();
        begin = new uint64_t(wordCount);
        end = begin + wordCount;
        copy(that.begin, that.end, begin);
        return *this;
    }

    // We also provide a move constructor, which makes it faster
    // to pass or return a BitSet by value.
    BitSet(BitSet&& that)
    {
        begin = that.begin;
        end = that.end;
        that.begin = 0;
        that.end = 0;
    }

    bool operator<(const BitSetPointer& that) const
    {
        return (*this) < that;
    }

};





// Count the number of mismatching bits between two bit vectors.
inline uint64_t ChanZuckerberg::ExpressionMatrix2::countMismatches(
    const BitSetPointer& x,
    const BitSetPointer& y)
{
    const uint64_t wordCount = x.wordCount();
    uint64_t mismatchCount = 0;
    for(uint64_t i = 0; i < wordCount; i++) {
        mismatchCount += __builtin_popcountll(x.begin[i] ^ y.begin[i]);
    }
    return mismatchCount;
}


#endif
