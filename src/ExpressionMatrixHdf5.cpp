// Functionality to read expression matrix data from an hdf5 file created
// by the 10X Genomics pipeline.

#include "ExpressionMatrix.hpp"
#include "hdf5.hpp"
#include "iterator.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;


/*******************************************************************************

Add cells from an hdf5 file in the format typically created by the 10X Genomics
pipeline.

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

The gene names are in data set gene_names. Data set genes is not used. It contains alternative
gene ids.

Cell names (or rather the corresponding barcodes) are in data set barcodes.

Data set shape is also not used.

*******************************************************************************/

void ExpressionMatrix::addCellsFromHdf5(
	const string& fileName,
	size_t maxTermCountForApproximateSimilarityComputation)
{

	try {
		// Open the file.
		const H5::H5File file(fileName, H5F_ACC_RDONLY);

		// Find the names of the top level groups in the file.
		// We expect only a single group.
		vector<string> groupNames;
		hdf5::getTopLevelGroups(file, groupNames);
		if(groupNames.size() != 1) {
			throw runtime_error("HDF5 file " + fileName + " should have one and only one group.");
		}
		string& groupName = groupNames.front();
		groupName = "/" + groupName + "/";

		// Read the gene names.
		const H5::DataSet geneNamesDataSet = file.openDataSet(groupName + "gene_names");
		vector<string> hdf5GeneNames;
		hdf5::read(geneNamesDataSet, hdf5GeneNames);
		cout << "Found " << hdf5GeneNames.size() << " genes." << endl;

		// Read the cell names.
		const H5::DataSet cellNamesDataSet = file.openDataSet(groupName + "barcodes");
		vector<string> hdf5CellNames;
		hdf5::read(cellNamesDataSet, hdf5CellNames);
		cout << "Found " << hdf5CellNames.size() << " cells." << endl;

	}
	catch (H5::Exception& e) {
		cout << "An error occurred while reading HDF5 file " << fileName << endl;
		cout << e.getDetailMsg() << endl;
		throw;
	}
	catch (std::exception& e) {
		cout << "An error occurred while reading HDF5 file " << fileName << endl;
		cout << e.what() << endl;
		throw;
	}

}
