#ifndef CZI_EXPRESSION_MATRIX2_BIT_SET_HPP
#define CZI_EXPRESSION_MATRIX2_BIT_SET_HPP


// A bare bone bitset class with just the functionality needed for operations
// on LSH signatures. Similar to boost::dynamic_bitset.

#include "CZI_ASSERT.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class BitSet;
        uint64_t countMismatches(const BitSet&, const BitSet&);
    }
}



class ChanZuckerberg::ExpressionMatrix2::BitSet {
public:

	BitSet(uint64_t bitCount)
	{
		CZI_ASSERT(bitCount > 0);
		data.resize(((bitCount-1ULL) >> 6ULL) + 1ULL, 0ULL);
	}



	// Get the bit at a given position.
	bool get(uint64_t bitPosition) const
	{
		// Find the word containing the bit.
		const uint64_t wordIndex = bitPosition >> 6ULL;
		const uint64_t& word = data[wordIndex];

		// Find the position of this bit in the word.
		const uint64_t bitPositionInWord = bitPosition & 63ULL;

		// Test the bit.
		const uint64_t mask = 1ULL << bitPositionInWord;
		return (word & mask) != 0ULL;

	}

	// Set a bit at a given position.
	void set(uint64_t bitPosition)
	{
		// Find the word containing the bit.
		const uint64_t wordIndex = bitPosition >> 6ULL;
		uint64_t& word = data[wordIndex];

		// Find the position of this bit in the word.
		const uint64_t bitPositionInWord = bitPosition & 63ULL;

		// Set the bit.
		word |= 1ULL << bitPositionInWord;

	}

	// Get a uint64_t containing the bits in a specified bit range.
	// This range must be entirely contained in one word (it cannot cross word boundaries).
	uint64_t getBits(uint64_t firstBitPosition, uint64_t bitCount) const
	{
		// Find the word containing the bit.
		const uint64_t firstBitWordIndex = firstBitPosition >> 6ULL;
		const uint64_t& word = data[firstBitWordIndex];

		// Find the position of the first bit in the word.
		const uint64_t firstBitPositionInWord = firstBitPosition & 63ULL;

		// Check that the entire range is contained in this word.
		CZI_ASSERT(firstBitPositionInWord + bitCount <= 64ULL);

		// Extract the bits.
		const uint64_t bitMask = (1ULL << bitCount) - 1ULL;
		return (word >> firstBitPositionInWord) & bitMask;
	}

	vector<uint64_t> data;


};


// Count the number of mismatching bits between two bit vectors.
// The two bit vectors should have the same length, but we don't check this for performance.
inline uint64_t ChanZuckerberg::ExpressionMatrix2::countMismatches(const BitSet& x, const BitSet& y)
{
	uint64_t mismatchCount = 0;
	const uint64_t blockCount = x.data.size();
	for(uint64_t i=0; i<blockCount; i++) {
		mismatchCount += __builtin_popcountll(x.data[i] ^ y.data[i]);
	}
	return mismatchCount;
}


#endif
