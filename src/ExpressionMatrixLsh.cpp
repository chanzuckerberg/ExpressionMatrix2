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
#include "lapack.hpp"
#include "nextPowerOfTwo.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/math/constants/constants.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <cmath>
#include "fstream.hpp"


// Generate the random unit LSH vectors.
// This uses the Marsaglia method, which consists of generating
// a vector with normally distributed components
// with zero mean and unit variance, then normalizing it, described in
// Marsaglia, G. "Choosing a Point from the Surface of a Sphere." Ann. Math. Stat. 43, 645-646, 1972.
// See http://mathworld.wolfram.com/HyperspherePointPicking.html
void ExpressionMatrix::generateLshVectors(
	size_t lshBandCount,
	size_t lshRowCount,
	unsigned int seed,
	vector< vector< vector<double> > >& lshVectors	// Indexed by [band][row][geneId]
    ) const
{
	// Prepare to generate normally vector distributed components.
	using RandomSource = boost::mt19937;
	using NormalDistribution = boost::normal_distribution<>;
	RandomSource randomSource(seed);
	NormalDistribution normalDistribution;
	boost::variate_generator<RandomSource, NormalDistribution> normalGenerator(randomSource, normalDistribution);

	// Allocate space for the LSH vectors.
	lshVectors.resize(lshBandCount, vector< vector<double> >(lshRowCount, vector<double>(geneCount())));



	// Triple loop over bands, rows, components (genes).
	for(size_t band=0; band<lshBandCount; band++) {
		for(size_t row=0; row<lshRowCount; row++) {
			for(GeneId geneId=0; geneId<geneCount(); geneId++) {
				lshVectors[band][row][geneId] = normalGenerator();
			}

			// Normalize the vector for this band and row.
			double sum2 = 0.;
			for(const double& x: lshVectors[band][row]) {
				sum2 += x*x;
			}
			const double factor = 1. / sqrt(sum2);
			for(double& x: lshVectors[band][row]) {
				x *= factor;
			}
		}
	}
}



// Approximate computation of the similarity between two cells using
// Locality Sensitive Hashing (LSH).
// Not to be used for code where performance is important,
// because it recompiutes the LSH vector every time.
double ExpressionMatrix::computeApproximateLshCellSimilarity(
	size_t lshBandCount,
	size_t lshRowCount,
	unsigned int seed,
	CellId cellId0,
	CellId cellId1) const
{

	// Generate LSH vectors.
	vector< vector< vector<double> > > lshVectors;
	generateLshVectors(lshBandCount, lshRowCount, seed, lshVectors);

	// Compute scalar products of each LSH vector with the expression counts of each cell.
	vector<double> scalarProducts0, scalarProducts1;
	for(size_t band=0; band<lshBandCount; band++) {
		for(size_t row=0; row<lshRowCount; row++) {
			scalarProducts0.push_back(computeExpressionCountScalarProduct(cellId0, lshVectors[band][row]));
			scalarProducts1.push_back(computeExpressionCountScalarProduct(cellId1, lshVectors[band][row]));
		}
	}

	// Count the number of scalar products with the same sign.
	size_t sameSignCount = 0;
	for(size_t i=0; i<scalarProducts0.size(); i++) {
		const double sign0 = std::copysign(1., scalarProducts0[i]);
		const double sign1 = std::copysign(1., scalarProducts1[i]);
		if(sign0 == sign1) {
			++sameSignCount;
		}
	}

	// Compute the estimated angle.
	const double angle =  boost::math::double_constants::pi * (1.-double(sameSignCount) / double(lshBandCount*lshRowCount));
	return std::cos(angle);
}



// Approximate computation of the angle between the expression vectors of two cells
// using Locality Sensitive Hashing (LSH).
// The approximate similarity can be computed as the cosine of this angle.
// This recomputes every time the scalar product of the cell expression vector
// with the LSH vectors.
double ExpressionMatrix::computeApproximateLshCellAngle(
	const vector< vector< vector<double> > >& lshVectors,
	CellId cellId0,
	CellId cellId1) const
{
	// Compute scalar products of each LSH vector with the expression counts of each cell.
	vector< pair<double, double> > scalarProducts;
	for(const auto& bandLshVectors: lshVectors) {
		for(const auto& lshVector: bandLshVectors) {
			scalarProducts.push_back(make_pair(
				computeExpressionCountScalarProduct(cellId0, lshVector),
				computeExpressionCountScalarProduct(cellId1, lshVector)));
		}
	}

	// Count the number of scalar products with the same sign.
	size_t sameSignCount = 0;
	for(size_t i=0; i<scalarProducts.size(); i++) {
		const double sign0 = std::copysign(1., scalarProducts[i].first);
		const double sign1 = std::copysign(1., scalarProducts[i].second);
		if(sign0 == sign1) {
			++sameSignCount;
		}
	}

	// Compute the estimated angle.
	const size_t totalCount = scalarProducts.size();
	const double angle =  boost::math::double_constants::pi * (1.-double(sameSignCount) / double(totalCount));
	return angle;
}



// Approximate computation of the angle between the expression vectors of two cells
// using Locality Sensitive Hashing (LSH).
double ExpressionMatrix::computeApproximateLshCellAngle(
	const BitSet& signature0,
	const BitSet& signature1,
	double bitCountInverse) const
{

	// Count the number of bits where the two signatures disagree.
	size_t oppositeSignCount = countMismatches(signature0, signature1);

	// Compute the estimated angle.
	const double angle =  boost::math::double_constants::pi * double(oppositeSignCount) * bitCountInverse;
	return angle;
}



// Compute the scalar product of an LSH vector with the normalized expression counts of a cell.
double ExpressionMatrix::computeExpressionCountScalarProduct(CellId cellId, const vector<double>& v) const
{
	// Sanity check on the vector being passed in.
	CZI_ASSERT(v.size() == geneCount());

	const Cell& cell = cells[cellId];
	const double mean = cell.sum1 / double(geneCount());
	const double sigmaInverse = 1./sqrt(cell.sum2/geneCount() - mean*mean);

	// Initialize the scalar product to what it would be
	// if all counts were zero.
	const double zeroContribution =  - mean * sigmaInverse;
	double scalarProduct = 0.;
	for(GeneId geneId=0; geneId<geneCount(); geneId++) {
		scalarProduct += v[geneId] * zeroContribution;
	}

	// Add the contribution of the non-zero expression counts for this cell.
	for(const auto& p: cellExpressionCounts[cellId]) {
		const GeneId geneId = p.first;
		const float& count = p.second;

		// Accumulate into the scalar product.
		scalarProduct += sigmaInverse * double(count) * v[geneId];
	}
	return scalarProduct;
}



// Write a csv file containing, for every pair of cells,
// the exact similarity and the similarity computed using LSH.
void ExpressionMatrix::writeLshSimilarityComparisonSlow(
	size_t lshBandCount,
	size_t lshRowCount,
	unsigned int seed
	) const
{
	// Generate LSH vectors.
	vector< vector< vector<double> > > lshVectors;
	generateLshVectors(lshBandCount, lshRowCount, seed, lshVectors);

	// Open the csv file
	ofstream csvOut("LshSimilarityComparison.csv");
	csvOut << "Cell0,Cell1,Exact,LSH,Delta\n";


    // Loop over all pairs.
    for(CellId cellId0=0; cellId0!=cellCount()-1; cellId0++) {
        if((cellId0%1) == 0) {
            cout << timestamp << "Working on cell " << cellId0 << " of " << cells.size() << endl;
        }
        for(CellId cellId1=cellId0+1; cellId1!=cellCount(); cellId1++) {
            const double exactSimilarity = computeCellSimilarity(cellId0, cellId1);
            const double approximateAngle = computeApproximateLshCellAngle(lshVectors, cellId0, cellId1);
            const double approximateSimilarity = std::cos(approximateAngle);
            const double delta = approximateSimilarity - exactSimilarity;
            csvOut << cellId0 << ",";
            csvOut << cellId1 << ",";
            csvOut << exactSimilarity << ",";
            csvOut << approximateSimilarity << ",";
            csvOut << delta << "\n";
            csvOut << flush;
        }
    }
}
void ExpressionMatrix::writeLshSimilarityComparison(
	size_t lshBandCount,
	size_t lshRowCount,
	unsigned int seed
	) const
{
	// Generate LSH vectors.
	vector< vector< vector<double> > > lshVectors;
	cout << timestamp << "Generating the LSH vectors." << endl;
	generateLshVectors(lshBandCount, lshRowCount, seed, lshVectors);

	// Orthogonalize the LSH vectors.
    // This did not seem to give any benefit, so I turned it off to eliminate the
    // dependency on Lapack.
	// cout << timestamp << "Orthogonalizing the LSH vectors." << endl;
	// orthogonalizeLshVectors(lshVectors, lshBandCount*lshRowCount);

	// Compute the LSH signatures of all cells.
	cout << timestamp << "Computing LSH signatures for all cells." << endl;
	vector<BitSet> signatures;
	computeCellLshSignatures(lshVectors, signatures);

	// Open the csv file
	ofstream csvOut("LshSimilarityComparison.csv");
	csvOut << "Cell0,Cell1,Exact,LSH,Delta\n";


    // Loop over all pairs.
	cout << timestamp << "Computing similarities for all cell pairs." << endl;
	const double bitCountInverse = 1. / double(lshBandCount * lshRowCount);
    for(CellId cellId0=0; cellId0!=cellCount()-1; cellId0++) {
        if((cellId0%100) == 0) {
            cout << timestamp << "Working on cell " << cellId0 << " of " << cells.size() << endl;
        }
        for(CellId cellId1=cellId0+1; cellId1!=cellCount(); cellId1++) {
            const double exactSimilarity = computeCellSimilarity(cellId0, cellId1);
            const double approximateAngle = computeApproximateLshCellAngle(signatures[cellId0], signatures[cellId1], bitCountInverse);
            const double approximateSimilarity = std::cos(approximateAngle);
            const double delta = approximateSimilarity - exactSimilarity;
            csvOut << cellId0 << ",";
            csvOut << cellId1 << ",";
            csvOut << exactSimilarity << ",";
            csvOut << approximateSimilarity << ",";
            csvOut << delta << "\n";
            csvOut << flush;
        }
    }
    cout << timestamp << endl;
}




// Given LSH vectors, compute the LSH signature of all cells.
// The LSH signature of a cell is a bit vector with one bit for each of the LSH vectors.
// Each bit is 1 if the scalar product of the cell expression vector
// (normalized to zero mean and unit variance) with the
// the LSH vector is positive, and 0 otherwise.
void ExpressionMatrix::computeCellLshSignatures(
	const vector< vector< vector<double> > >& lshVectors,
	vector<BitSet>& signatures
	) const
{
	// Count the total number of lsh vectors, assuming that all the bands have the same
	// number of rows.
	CZI_ASSERT(!lshVectors.empty());
	const size_t lshBandCount = lshVectors.size();
	CZI_ASSERT(!lshVectors.front().empty());
	const size_t lshRowCount = lshVectors.front().size();
	const size_t bitCount = lshBandCount * lshRowCount;

	// Compute the sum of the components of each lsh vector.
	// It is needed below to compute the contribution of the
	// expression counts that are zero.
	vector<double> lshVectorsSums(bitCount);
	size_t index = 0;
	for(size_t band=0; band<lshBandCount; band++) {
		const auto& bandVectors = lshVectors[band];
		CZI_ASSERT(bandVectors.size() == lshRowCount);	// All bands must have the same number of rows.
		for(size_t row=0; row<lshRowCount; row++, index++) {
			const auto& rowVector = bandVectors[row];
			CZI_ASSERT(rowVector.size() == geneCount());
			lshVectorsSums[index] = std::accumulate(rowVector.begin(), rowVector.end(), 0.);
		}
	}

	// Initialize the signatures to all zero bits.
	signatures.resize(cellCount(), BitSet(bitCount));

	// Vector to hold the scalar products of the normalized expression vector of a cell
	// with each of the LSH vectors.
	vector<double> scalarProducts(bitCount);

	// Loop over all cells.
	for(CellId cellId=0; cellId<cellCount(); cellId++) {
		if((cellId>0) && ((cellId % 100)==0)) {
			cout << timestamp << "Working on cell " << cellId << endl;
		}

		// Compute the mean and standard deviation for this cell.
		const Cell& cell = cells[cellId];
		const double mean = cell.sum1 / double(geneCount());
		const double sigmaInverse = 1./sqrt(cell.sum2/geneCount() - mean*mean);

		// Initialize the scalar products to what they would be if all counts were zero.
		const double zeroContribution =  - mean * sigmaInverse;
		for(size_t index=0; index<bitCount; index++) {
			scalarProducts[index] = lshVectorsSums[index] * zeroContribution;
		}

		// Add the contributions of the non-zero expression counts for this cell.
		for(const auto& p: cellExpressionCounts[cellId]) {
			const GeneId geneId = p.first;
			const float& count = p.second;
			const double scaledCount = sigmaInverse * double(count);

			// Accumulate into the scalar products.
			size_t index = 0;
			for(size_t band=0; band<lshBandCount; band++) {
				const auto& bandVectors = lshVectors[band];
				for(size_t row=0; row<lshRowCount; row++, index++) {
					CZI_ASSERT(index < scalarProducts.size());
					scalarProducts[index] += scaledCount * bandVectors[row][geneId];
				}
			}
		}

		// Set to 1 the signature bits corresponding to positive scalar products.
		auto& cellSignature = signatures[cellId];
		for(size_t index=0; index<bitCount; index++) {
			if(scalarProducts[index] >0) {
				cellSignature.set(index);
			}
		}
	}
}




// Same as above, but only for a set of cells given in a vector of cell ids (cell set).
void ExpressionMatrix::computeCellLshSignatures(
	const vector< vector< vector<double> > >& lshVectors,
	const MemoryMapped::Vector<CellId>& cellSet,
	vector<BitSet>& signatures
	) const
{
	// Count the total number of lsh vectors, assuming that all the bands have the same
	// number of rows.
	CZI_ASSERT(!lshVectors.empty());
	const size_t lshBandCount = lshVectors.size();
	CZI_ASSERT(!lshVectors.front().empty());
	const size_t lshRowCount = lshVectors.front().size();
	const size_t bitCount = lshBandCount * lshRowCount;

	// Compute the sum of the components of each lsh vector.
	// It is needed below to compute the contribution of the
	// expression counts that are zero.
	vector<double> lshVectorsSums(bitCount);
	size_t index = 0;
	for(size_t band=0; band<lshBandCount; band++) {
		const auto& bandVectors = lshVectors[band];
		CZI_ASSERT(bandVectors.size() == lshRowCount);	// All bands must have the same number of rows.
		for(size_t row=0; row<lshRowCount; row++, index++) {
			const auto& rowVector = bandVectors[row];
			CZI_ASSERT(rowVector.size() == geneCount());
			lshVectorsSums[index] = std::accumulate(rowVector.begin(), rowVector.end(), 0.);
		}
	}

	// Initialize the signatures to all zero bits.
	signatures.resize(cellSet.size(), BitSet(bitCount));

	// Vector to hold the scalar products of the normalized expression vector of a cell
	// with each of the LSH vectors.
	vector<double> scalarProducts(bitCount);

	// Loop over all cells.
	for(CellId localCellId=0; localCellId<cellSet.size(); localCellId++) {
		if((localCellId>0) && ((localCellId % 100)==0)) {
			cout << timestamp << "Working on cell " << localCellId << " of " << cellSet.size() << endl;
		}

		// Compute the mean and standard deviation for this cell.
		const CellId globalCellId = cellSet[localCellId];
		const Cell& cell = cells[globalCellId];
		const double mean = cell.sum1 / double(geneCount());
		const double sigmaInverse = 1./sqrt(cell.sum2/geneCount() - mean*mean);

		// Initialize the scalar products to what they would be if all counts were zero.
		const double zeroContribution =  - mean * sigmaInverse;
		for(size_t index=0; index<bitCount; index++) {
			scalarProducts[index] = lshVectorsSums[index] * zeroContribution;
		}

		// Add the contributions of the non-zero expression counts for this cell.
		for(const auto& p: cellExpressionCounts[globalCellId]) {
			const GeneId geneId = p.first;
			const float& count = p.second;
			const double scaledCount = sigmaInverse * double(count);

			// Accumulate into the scalar products.
			size_t index = 0;
			for(size_t band=0; band<lshBandCount; band++) {
				const auto& bandVectors = lshVectors[band];
				for(size_t row=0; row<lshRowCount; row++, index++) {
					CZI_ASSERT(index < scalarProducts.size());
					scalarProducts[index] += scaledCount * bandVectors[row][geneId];
				}
			}
		}

		// Set to 1 the signature bits corresponding to positive scalar products.
		auto& cellSignature = signatures[localCellId];
		for(size_t index=0; index<bitCount; index++) {
			if(scalarProducts[index] >0) {
				cellSignature.set(index);
			}
		}
	}
}



// Write to a csv file statistics of the cell LSH signatures..
void ExpressionMatrix::writeLshSignatureStatistics(size_t bitCount, const vector<BitSet>& signatures) const
{
	ofstream csvOut("LshSignatureStatistics.csv");
	csvOut << "Bit,Set,Unset,Total\n";

	for(size_t i=0; i<bitCount; i++) {

		// Count the numebr of cells that have this bit set.
		size_t setCount = 0;
		for(CellId cellId=0; cellId<cellCount(); cellId++) {
			if(signatures[cellId].get(i)) {
				++setCount;
			}
		}
		const size_t unsetCount = cellCount() - setCount;

		csvOut << i << "," << setCount << "," << unsetCount << "," << cellCount() << "\n";

	}

}




// Find similar cell pairs by looping over all pairs
// and using an LSH approximation to compute the similarity between two cells.
// See the beginning of ExpressionMatrixLsh.cpp for more information.
// Like findSimilarPairs0, this is also O(N**2) slow. However
// the coefficient of the N**2 term is much lower (around 15 ns/pair when nothing gets stored), at a cost of
// additional O(N) work (typically 30 ms per cell for lshCount=1024).
// As a result, this can be much faster for large numbers of cells.
// The error of the approximation is controlled by lshCount.
// The maximum standard deviation of the computed similarity is (pi/2)/sqrt(lshCount),
// or about 0.05 for lshCount=1024.
// The standard deviation decreases as the similarity increases. It becomes
// zero when the similarity is 1. For similarity 0.5, the standard deviation is 82%
// of the standard deviation at similarity 0.
void ExpressionMatrix::findSimilarPairs1(
	const string& cellSetName,	// The name of the cell set to be used.
    const string& name,         // The name of the SimilarPairs object to be created.
    size_t k,                   // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold, // The minimum similarity for a pair to be stored.
	size_t lshCount,			// The number of LSH functions (hyperplanes) to be used.
	unsigned int seed 			// The seed used to generate the random hyperplanes.
	)
{
	// Sanity check.
	CZI_ASSERT(similarityThreshold <= 1.);

	// Locate the cell set.
	const auto& it = cellSets.cellSets.find(cellSetName);
	if(it == cellSets.cellSets.end()) {
		throw runtime_error("Cell set " + cellSetName + " does not exist.");
	}
	const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
	if(cellSet.size() == 0) {
		cout << "Cell set " << cellSetName << " is empty. Skipping findSimilarPairs1." << endl;
		return;
	}

	// Generate LSH vectors.
	const size_t lshBandCount = lshCount;
	const size_t lshRowCount = 1;
	vector< vector< vector<double> > > lshVectors;
	generateLshVectors(lshBandCount, lshRowCount, seed, lshVectors);

	// Compute the LSH signatures of all cells.
	cout << timestamp << "Computing LSH signatures for all cells." << endl;
	vector<BitSet> signatures;
	computeCellLshSignatures(lshVectors, cellSet, signatures);

	// Write signature statistics to a csv file.
	// writeLshSignatureStatistics(lshCount, signatures);

    // Create the SimilarPairs object where we will store the pairs.
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + name, k, cellSet);

    // Compute the angle threshold (in radians) corresponding to this similarity threshold.
    const double angleThreshold = std::acos(similarityThreshold);

    // Compute the corresponding threshold on the number of mismatching
    // signature bits.
    const size_t mismatchCountThreshold = size_t(double(lshCount) * angleThreshold / boost::math::double_constants::pi);

    // Compute the similarity (cosine of the angle) corresponding to every number of mismatching bits
    // up to the threshold.
    vector<double> similarityTable(mismatchCountThreshold + 1);
    for(size_t mismatchingBitCount=0; mismatchingBitCount<=mismatchCountThreshold; mismatchingBitCount++) {
    	const double angle = double(mismatchingBitCount) * boost::math::double_constants::pi / double(lshCount);
    	CZI_ASSERT(mismatchingBitCount < similarityTable.size());
    	similarityTable[mismatchingBitCount] = std::cos(angle);
    }



    // Loop over all pairs. This is much faster than findSimilarPairs0
    // (around 15 ns per pair when using LSH vectors of 1024 bits), but
    // still scales like the square of the number of cells in the cell set.
	cout << timestamp << "Begin computing similarities for all cell pairs." << endl;
    for(CellId localCellId0=0; localCellId0!=cellSet.size()-1; localCellId0++) {
        if(localCellId0>0 && ((localCellId0%10000) == 0)) {
            cout << timestamp << "Working on cell " << localCellId0 << " of " << cellSet.size() << endl;
        }
        const BitSet& signature0 = signatures[localCellId0];
        for(CellId localCellId1=localCellId0+1; localCellId1!=cellSet.size(); localCellId1++) {
            const BitSet& signature1 = signatures[localCellId1];

        	// Count the number of bits where the signatures of these two cells disagree.
        	size_t mismatchingBitCount = countMismatches(signature0, signature1);

            // If the similarity is sufficient, pass it to the SimilarPairs container,
            // which will make the decision whether to store it, depending on the
            // number of pairs already stored for cellId0 and cellId1.
            if(mismatchingBitCount <= mismatchCountThreshold) {
            	CZI_ASSERT(mismatchingBitCount < similarityTable.size());
                similarPairs.add(localCellId0, localCellId1, similarityTable[mismatchingBitCount]);
            }
        }
    }



    // Sort the similar pairs for each cell by decreasing similarity.
	cout << timestamp << "Sorting pairs." << endl;
    similarPairs.sort();
	cout << timestamp << "Done sorting pairs." << endl;

}



// Find similar cell pairs using LSH, without looping over all pairs.
// See the beginning of ExpressionMatrixLsh.cpp for more information.
// This implementation requires lshRowCount to be a power of 2 not greater than 64.
void ExpressionMatrix::findSimilarPairs2(
	const string& cellSetName,	// The name of the cell set to be used.
    const string& name,         // The name of the SimilarPairs object to be created.
    size_t k,                   // The maximum number of similar pairs to be stored for each cell.
    double similarityThreshold, // The minimum similarity for a pair to be stored.
	size_t lshBandCount,		// The number of LSH bands, each generated using lshRowCount LSH vectors.
	size_t lshRowCount,         // The number of LSH vectors in each of the lshBandCount LSH bands.
	unsigned int seed,			// The seed used to generate the LSH vectors.
	double loadFactor           // Of the hash table used to assign cells to bucket.
	)
{
	// Sanity check.
	CZI_ASSERT(similarityThreshold <= 1.);

	// Locate the cell set.
	const auto& it = cellSets.cellSets.find(cellSetName);
	if(it == cellSets.cellSets.end()) {
		throw runtime_error("Cell set " + cellSetName + " does not exist.");
	}
	const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
	if(cellSet.size() == 0) {
		cout << "Cell set " << cellSetName << " is empty. Skipping findSimilarPairs1." << endl;
		return;
	}

    // Compute the angle threshold (in radians) corresponding to this similarity threshold.
    const double angleThreshold = std::acos(similarityThreshold);

    // Compute the corresponding threshold on the number of mismatching
    // signature bits.
    const size_t lshCount = lshBandCount * lshRowCount;
    const size_t mismatchCountThreshold = size_t(double(lshCount) * angleThreshold / boost::math::double_constants::pi);

    // Compute the similarity (cosine of the angle) corresponding to every number of mismatching bits
    // up to the threshold.
    vector<double> similarityTable(mismatchCountThreshold + 1);
    for(size_t mismatchingBitCount=0; mismatchingBitCount<=mismatchCountThreshold; mismatchingBitCount++) {
    	const double angle = double(mismatchingBitCount) * boost::math::double_constants::pi / double(lshCount);
    	CZI_ASSERT(mismatchingBitCount < similarityTable.size());
    	similarityTable[mismatchingBitCount] = std::cos(angle);
    }


	// Generate LSH vectors.
	vector< vector< vector<double> > > lshVectors;
	generateLshVectors(lshBandCount, lshRowCount, seed, lshVectors);

	// Compute the LSH signatures of all cells.
	cout << timestamp << "Computing LSH signatures for all cells." << endl;
	vector<BitSet> signatures;
	computeCellLshSignatures(lshVectors, cellSet, signatures);

    // Create the SimilarPairs object where we will store the pairs.
    SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + name, k, cellSet);

    // Vector of vectors used to contain the cells assigned to each bucket.
    const uint64_t bucketCount = nextPowerOfTwoGreaterThanOrEqual(uint64_t(double(cellCount()) / loadFactor));
    const uint64_t bucketMask = bucketCount - 1ULL;
    MemoryMapped::VectorOfVectors<CellId, uint64_t> buckets;
    buckets.createNew(directoryName + "/SimilarPairsBuckets-" + name);



    // Loop over the LSH bands.
	uint64_t checkedPairCount = 0ULL;
    for(size_t band=0; band<lshBandCount; band++) {
    	cout << timestamp << "Begin working on band " << band << " of " << lshBandCount << endl;

    	// The index of the first bit of this band.
    	const uint64_t firstBit = band * lshRowCount;

    	// Assign cells to buckets. In pass 1 we just count the number of cells
    	// in each bucket.
    	buckets.beginPass1(bucketCount);
    	for(CellId cellId=0; cellId<cellCount(); cellId++) {
    		const BitSet& signature = signatures[cellId];
    		const uint64_t signatureSlice = signature.getBits(firstBit, lshRowCount);
    		const uint64_t bucket = MurmurHash64A(&signatureSlice, sizeof(signatureSlice), 231) & bucketMask;
    		buckets.incrementCount(bucket);
    	}

    	// Assign cells to buckets. In pass 2 we actually store the cells in buckets.
    	buckets.beginPass2();
    	for(CellId cellId=0; cellId<cellCount(); cellId++) {
    		const BitSet& signature = signatures[cellId];
    		const uint64_t signatureSlice = signature.getBits(firstBit, lshRowCount);
    		const uint64_t bucket = MurmurHash64A(&signatureSlice, sizeof(signatureSlice), 231) & bucketMask;
    		buckets.store(bucket, cellId);
    	}
    	buckets.endPass2();



    	// Now loop over buckets.
    	// For each bucket, check all pairs of cells in the bucket.
    	uint64_t maxBucketSize = 0ULL;
    	uint64_t bandCheckedPairCount = 0ULL;
    	for(uint64_t bucketIndex=0; bucketIndex<bucketCount; bucketIndex++) {
    		const auto& bucket = buckets[bucketIndex];
    		const auto bucketSize = bucket.size();
    		maxBucketSize = max(maxBucketSize, bucketSize);

    		// If the bucket contains less than two cells, there is nothing to do.
    		if(bucketSize < 2) {
    			continue;
    		}

    		// Check all pairs in this bucket.
    		for(size_t i0=1; i0<bucket.size(); i0++) {
    			const CellId cellId0 = bucket[i0];
    			const auto& signature0 = signatures[cellId0];
    			for(size_t i1=0; i1<i0; i1++, ++bandCheckedPairCount) {
        			const CellId cellId1 = bucket[i1];
        			const auto& signature1 = signatures[cellId1];

        			// Count the number of bits where the signatures of these two cells disagree.
		        	size_t mismatchingBitCount = countMismatches(signature0, signature1);

		            // If the similarity is sufficient, pass it to the SimilarPairs container,
		            // which will make the decision whether to store it, depending on the
		            // number of pairs already stored for cellId0 and cellId1.
		            if(mismatchingBitCount <= mismatchCountThreshold) {
		            	CZI_ASSERT(mismatchingBitCount < similarityTable.size());
		                similarPairs.add(cellId0, cellId1, similarityTable[mismatchingBitCount]);
		            }
    			}
    		}
    	}
    	cout << "Band " << band << ": maximum bucket size " << maxBucketSize;
    	cout << ". Checked " << bandCheckedPairCount << " pairs." << endl;
    	checkedPairCount += bandCheckedPairCount;

    }
    cout << "Total number of pairs checked is " << checkedPairCount << endl;


    // Sort the similar pairs for each cell by decreasing similarity.
	cout << timestamp << "Sorting pairs." << endl;
    similarPairs.sort();
	cout << timestamp << "Done sorting pairs." << endl;
}



#if 0
// Orthogonalize the LSH vectors in groups of k.
// Ji et al. (2012) have shown that orthogonalization of
// groups of LSH vectors can result in significant reduction of the
// variance of the distance estimate provided by LSH.
// See J. Ji, J. Li, S. Yan, B. Zhang, and Q. Tian,
// Super-Bit Locality-Sensitive Hashing, In NIPS, pages 108–116, 2012.
// https://pdfs.semanticscholar.org/64d8/3ccbcb1d87bfafee57f0c2d49043ee3f565b.pdf
// This did not seem to give any benefit, so I turned it off to eliminate the
// dependency on Lapack.
void ExpressionMatrix::orthogonalizeLshVectors(
	vector< vector< vector<double> > >& lshVectors,
	size_t k
    ) const
{

	// Create a vector of pair(band, row) for all the LSH vectors,
	// We will then partition that in groups of k.
	const size_t lshBandCount = lshVectors.size();
	CZI_ASSERT(lshBandCount > 0);
	const size_t lshRowCount = lshVectors.front().size();
	CZI_ASSERT(lshRowCount > 0);
	vector< pair<size_t, size_t> > lshVectorsInfo;
	for(size_t band=0; band<lshBandCount; band++) {
		CZI_ASSERT(lshVectors[band].size() == lshRowCount);
		for(size_t row=0; row<lshRowCount; row++) {
			CZI_ASSERT(lshVectors[band][row].size() == geneCount());
			lshVectorsInfo.push_back(make_pair(band, row));
		}
	}



	// Matrix to hold a group of k LSH vectors to be orthogonalized.
	// Each vector is a column of the matrix, of length geneCount(),
	// stored contiguously in memory. This is the format required
	// by Lapack routine dgeqrf.
	vector<double> a(geneCount() * k);

	// Work vectors for Lapack calls.
	vector<double> tau(k);
	const size_t lWork = k * 128;
	vector<double> work(lWork);



	// Loop over groups of k LSH vectors.
	// Orthogonalize each group.
	for(size_t begin=0; begin<lshVectorsInfo.size(); begin+=k) {
		const size_t end = min(begin+k, lshVectorsInfo.size());

		// The number of vectors we are orthogonalizing in this iteration.
		// Always equal to k, except possibly for the last group.
		const size_t n = end - begin;

		// Copy the LSH vectors into the columns of matrix a.
		size_t aIndex = 0ULL;
		for(size_t i=begin; i!=end; i++) {
			const auto& p = lshVectorsInfo[i];
			const size_t band = p.first;
			const size_t row = p.second;
			const vector<double>& v = lshVectors[band][row];
			CZI_ASSERT(v.size() == geneCount());
			copy(v.begin(), v.end(), a.begin()+aIndex);
			aIndex += geneCount();
		}
		CZI_ASSERT(aIndex == n * geneCount());



		// Use Lapack to orthogonalize the vectors stored in a.
		// This uses a Householder QR factorization, which is
		// numerically much more stable than the simple
		// Gram-Schmidt orthogonalization procedure.
	    int info = 0;
	    dgeqrf_(int(geneCount()), int(n), &a.front(), int(geneCount()), &tau.front(), &work.front(), int(lWork), info);
	    CZI_ASSERT(info == 0);
	    dorgqr_(int(geneCount()), int(n), int(n), &a.front(), int(geneCount()), &tau.front(), &work.front(), int(lWork), info);
	    CZI_ASSERT(info == 0);


		// Copy back the columns of matrix a,
		// which now contain the orthogonalized LSH vectors.
		aIndex = 0ULL;
		for(size_t i=begin; i!=end; i++) {
			const auto& p = lshVectorsInfo[i];
			const size_t band = p.first;
			const size_t row = p.second;
			vector<double>& v = lshVectors[band][row];
			CZI_ASSERT(v.size() == geneCount());
			copy(a.begin()+aIndex, a.begin()+aIndex+geneCount(), v.begin());
			aIndex += geneCount();
		}
		CZI_ASSERT(aIndex == n * geneCount());
	}
}
#endif

