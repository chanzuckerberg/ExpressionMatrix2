#ifndef CZI_EXPRESSION_MATRIX2_SKIP_HDF5
// Functions used to read hdf5 files containing expression matrices
// created by the 10X Genomics pipeline.

// See hdf5.hpp for more information.

#include "hdf5.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;
using namespace hdf5;

#include "iostream.hpp"



// Read a vector of strings from an hdf5 data set with rank 1.
void hdf5::read(
	const H5::DataSet& dataSet,
    vector<string>& data)
{
    // Check that the data set contains strings.
    if(dataSet.getTypeClass() != H5T_STRING) {
        throw runtime_error("Unexpected type class while reading strings from hdf5 data set.");
    }

    // Get the data space.
    const DataSpace dataSpace = dataSet.getSpace();

    // Check that the data space has rank 1 (contains a vector).
    if(dataSpace.getSimpleExtentNdims() != 1) {
        throw runtime_error("Unexpected rank while reading strings from hdf5 data set.");
    }

    // Get the number of strings.
    hsize_t n;
    dataSpace.getSimpleExtentDims(&n, 0);

    // Get the string type and its maximum length.
    const StrType strType(dataSet);
    const size_t stringMaximumLength = strType.getSize();

    // Now we know how many bytes to allocate and we can read the data in.
    vector<char> buffer(n * stringMaximumLength);
    dataSet.read(&buffer[0], strType);

    // Store the strings.
    data.resize(n);
    const char* bufferBegin = &buffer.front();
    for(size_t i=0; i<n; i++) {
        const char* address = bufferBegin + stringMaximumLength*i;
        string& s = data[i];

        // The string in the buffer is not null terminated
        // if it is of maximum length (Fortran style), so we have to be a bit careful.
        s = string(address, stringMaximumLength);
        const auto nullPosition = s.find(char(0));
        if(nullPosition != string::npos) {
            s.resize(nullPosition);
        }

    }

}



// Get into a vector of strings the names of the top level groups in the file.
// The helper function is called by H5Literate.
static herr_t getTopLevelGroupsHelper(
	hid_t,
	const char* name,
	const H5L_info_t*,
	void* p)
{
	vector<string>& groupNames = *static_cast< vector<string>* >(p);
    groupNames.push_back(name);
    return 0;
}
void hdf5::getTopLevelGroups(
	const H5::H5File& file,
	vector<string>& groupNames)
{
	H5Literate(file.getId(), H5_INDEX_NAME, H5_ITER_INC, 0, getTopLevelGroupsHelper, &groupNames);
}

#endif
