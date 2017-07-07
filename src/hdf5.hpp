#ifndef CZI_EXPRESSION_MATRIX2_HDF5_HPP
#define CZI_EXPRESSION_MATRIX2_HDF5_HPP



/*******************************************************************************

Functions used to read hdf5 files containing expression matrices
created by the 10X Genomics pipeline.

See here for information on the hdf5 C++ API:
https://support.hdfgroup.org/HDF5/doc/cpplus_RM/index.html

On Ubuntu 16.04, this requires the following packages:
libhdf5-10       (provides libhdf5_serial_hl.so)
libhdf5-cpp-11   (provides libhdf5_cpp.so)
libhdf5-dev      (only needed for building - provides the necessary include files)
hdf5-tools       (optional - command line tools to manipulate hdf5 files)

This requires adding the following directory to the include path for compilation:
/usr/include/hdf5/serial

It also requires the following libraries to be added to the link line:
hdf5_cpp
hdf5_serial

See here for an example of a typical hdf5 file that we want to process:
http://cf.10xgenomics.com/samples/cell-exp/1.3.0/t_3k_4k_aggregate/t_3k_4k_aggregate_filtered_gene_bc_matrices_h5.h5

You can use command hdf5ls (from package hdf5-tools) to list the contents of the file:

h5ls -r t_3k_4k_aggregate_filtered_gene_bc_matrices_h5.h5
/                        Group
/GRCh38                  Group
/GRCh38/barcodes         Dataset {8083}
/GRCh38/data             Dataset {8863417/Inf}
/GRCh38/gene_names       Dataset {33694}
/GRCh38/genes            Dataset {33694}
/GRCh38/indices          Dataset {8863417/Inf}
/GRCh38/indptr           Dataset {8084/8192}
/GRCh38/shape            Dataset {2/16384}

Data sets data, indices, and indptr contain the expression matrix in sparse format, by cell.

The gene names are in data set gene_names. Data set genes contains alternative
gene ids.

Cell names (or rather the corresponding barcodes) are in data set barcodes.

*******************************************************************************/

#include "stdexcept.hpp"
#include "string.hpp"
#include "vector.hpp"

// Include file for the hdf5 C++ API.
#include "H5Cpp.h"



namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
    namespace hdf5 {

    	// Add to this namespace the H5 namespace, which contains the hdf5 C++ API.
    	using namespace H5;

    	// Get into a vector of strings the names of the top level groups in the file.
    	void getTopLevelGroups(
			const H5::H5File&,
			vector<string>& groupNames);

		// Read a vector of integers from an hdf5 data set with rank 1.
		template<class Int> inline void read(
			const H5::DataSet&,
			vector<Int>& data);

		// Read a vector of integers from a portion of an hdf5 data set with rank 1.
		// The portion to read is specified by:
		// - offset: the index of the first element to read.
		// - n: the number of elements to read.
		template<class Int> inline void read(
			const H5::DataSet&,
			hsize_t offset,
			hsize_t n,
			vector<Int>& data);

		// Read a vector of strings from an hdf5 data set with rank 1.
		void read(
			const H5::DataSet&,
			vector<string>& data);

		}
	}
}



// Read a vector of integers from an hdf5 data set with rank 1.
template<class Int> inline void ChanZuckerberg::ExpressionMatrix2::hdf5::read(
	const H5::DataSet& dataSet,
    vector<Int>& data)
{
    // Sanity check on the integer type.
    static_assert(std::is_integral<Int>::value, "Instantiation of hd5::read requires an integer type.");

    // Check that the data set contains integers.
    if(dataSet.getTypeClass() != H5T_INTEGER) {
        throw runtime_error("Unexpected type class while reading integers from hdf5 data set.");
    }

    // Get the data space.
    const DataSpace dataSpace = dataSet.getSpace();

    // Check that the data space has rank 1 (contains a vector).
    if(dataSpace.getSimpleExtentNdims() != 1) {
        throw runtime_error("Unexpected rank while reading integers from hdf5 data set.");
    }

    // Get the number of integers.
    hsize_t n;
    dataSpace.getSimpleExtentDims(&n, 0);

    // Get the integer type and its length.
    const IntType intType(dataSet);
    if(intType.getSize() != sizeof(Int)) {
        throw runtime_error("Unexpected integer size while reading integers from hdf5 data set.");
    }

    // Now we know how many bytes to allocate and we can read the data in directly into the vector.
    data.resize(n);
    dataSet.read(&data.front(), intType);


}



// Read a vector of integers from a portion of an hdf5 data set with rank 1.
// The portion to read is specified by:
// - offset: the index of the first element to read.
// - n: the number of elements to read.
template<class Int> inline void ChanZuckerberg::ExpressionMatrix2::hdf5::read(
	const H5::DataSet& dataSet,
    hsize_t offset,
    hsize_t n,
    vector<Int>& data)
{
    // Sanity check on the integer type.
    static_assert(std::is_integral<Int>::value, "Instantiation of hd5::read requires an integer type.");

    // Check that the data set contains integers.
    if(dataSet.getTypeClass() != H5T_INTEGER) {
        throw runtime_error("Unexpected type class while reading integers from hdf5 data set.");
    }

    // Get the data space.
    DataSpace dataSpace = dataSet.getSpace();

    // Check that the data space has rank 1 (contains a vector).
    if(dataSpace.getSimpleExtentNdims() != 1) {
        throw runtime_error("Unexpected rank while reading integers from hdf5 data set.");
    }

    // Get the integer type and its length.
    const IntType intType(dataSet);
    if(intType.getSize() != sizeof(Int)) {
        throw runtime_error("Unexpected integer size while reading integers from hdf5 data set.");
    }

    // Define the hyperslab we want to read.
    dataSpace.selectHyperslab(H5S_SELECT_SET, &n, &offset);

    // Define the memory data space that will cause the entire buffer to be filled.
    DataSpace memoryDataSpace(1, &n);

    // Read the data in directly into the vector.
    data.resize(n);
    dataSet.read(&data.front(), intType, memoryDataSpace, dataSpace);
}



#endif
