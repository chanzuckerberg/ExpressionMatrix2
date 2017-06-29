#include "CZI_ASSERT.hpp"
#include "regressionCoefficient.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



// Compute the Pearson correlation coefficient of two vectors.
// See https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
double ChanZuckerberg::ExpressionMatrix2::regressionCoefficient(
	const vector<double>& X,
	const vector<double>& Y
	)
{
	// Get the number of points.
	const size_t n = X.size();

	// The two vectors must have the same length.
	CZI_ASSERT(Y.size() == n);

	// Compute the required sums.
	double sx = 0.;
	double sy = 0.;
	double sxx = 0.;
	double syy = 0.;
	double sxy = 0.;
	for(size_t i=0; i<n; i++) {
		const double x = X[i];
		const double y = Y[i];
		sx += x;
		sy += y;
		sxx += x*x;
		syy += y*y;
		sxy += x*y;
	}


	// Compute the correlation coefficient.
	const double nDouble = double(n);
	const double numerator = nDouble*sxy - sx*sy;
	const double denominator = sqrt((nDouble*sxx-sx*sx) * (nDouble*syy-sy*sy));
	return numerator / denominator;
}


