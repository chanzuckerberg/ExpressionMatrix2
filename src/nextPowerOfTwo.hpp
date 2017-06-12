#ifndef CZI_EXPRESSION_MATRIX2_ROUND_TO_NEXT_POWER_OF_TWO_HPP
#define CZI_EXPRESSION_MATRIX2_ROUND_TO_NEXT_POWER_OF_TWO_HPP

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        inline uint64_t nextPowerOfTwoGreaterThan(uint64_t x)
        {
            if(x == 0ULL) {
                return 1ULL;
            } else {
                const uint64_t leadingZeroBitCount = __builtin_clzl(x);
                return 1ULL << (64ULL - leadingZeroBitCount);
            }

        }

        inline uint64_t nextPowerOfTwoGreaterThanOrEqual(uint64_t x)
        {
            return nextPowerOfTwoGreaterThan(x - 1ULL);
        }

    }
}




#endif
