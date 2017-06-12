#ifndef CZI_EXPRESSION_MATRIX2_DEDUPLICATE_HPP
#define CZI_EXPRESSION_MATRIX2_DEDUPLICATE_HPP

#include "algorithm.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        template<class T> void deduplicate(vector<T>& v)
        {
            sort(v.begin(), v.end());
            v.resize(unique(v.begin(), v.end()) - v.begin());
        }
    }
}

#endif

