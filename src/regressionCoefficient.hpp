#ifndef CZI_EXPRESSION_MATRIX2_REGRESSION_COEFFICIENT_HPP
#define CZI_EXPRESSION_MATRIX2_REGRESSION_COEFFICIENT_HPP

#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {


    	// Compute the Pearson correlation coefficient of two vectors.
    	// See https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    	double regressionCoefficient(
    		const vector<double>& x,
			const vector<double>& y
			);

    }
}



#endif

