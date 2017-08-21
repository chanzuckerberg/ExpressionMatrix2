// Functionality to read expression matrix data from an hdf5 file created
// by the 10X Genomics pipeline.

#include "ExpressionMatrix.hpp"
#include "hdf5.hpp"
#include "iterator.hpp"
#include "timestamp.hpp"
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

The gene names are in data set genes. Data set gene_names cannot used because it can
contain duplications.

Cell names (or rather the corresponding barcodes) are in data set barcodes.

Data set shape is also not used.

*******************************************************************************/

void ExpressionMatrix::addCellsFromHdf5(const string& fileName)
{

    try {
        // Open the file.
        const H5::H5File file(fileName, H5F_ACC_RDONLY);

        // Find the names of the top level groups in the file.
        // We expect only a single group.
        vector<string> groupNames;
        hdf5::getTopLevelGroups(file, groupNames);
        if (groupNames.size() != 1) {
            throw runtime_error("HDF5 file " + fileName + " should have one and only one group.");
        }
        string& groupName = groupNames.front();
        groupName = "/" + groupName + "/";

        // Read the gene names.
        const H5::DataSet geneNamesDataSet = file.openDataSet(groupName + "genes");
        vector<string> hdf5GeneNames;
        hdf5::read(geneNamesDataSet, hdf5GeneNames);
        cout << "Found " << hdf5GeneNames.size() << " genes." << endl;

        // Check for duplications in the gene names.
        {
            vector<string> sortedHdf5GeneNames = hdf5GeneNames;
            sort(sortedHdf5GeneNames.begin(), sortedHdf5GeneNames.end());
            for (size_t i = 1; i < sortedHdf5GeneNames.size(); i++) {
                if (sortedHdf5GeneNames[i - 1] == sortedHdf5GeneNames[i]) {
                    const string& duplicateGeneName = sortedHdf5GeneNames[i - 1];
                    cout << "Duplicate gene name " + duplicateGeneName << " in hdf5 file " << fileName << endl;
                    cout << "This gene appears at the following positions in HDF5 data set " << groupName << "genes:";
                    for (size_t j = 0; j < hdf5GeneNames.size(); j++) {
                        if (hdf5GeneNames[j] == duplicateGeneName) {
                            cout << " " << j;
                        }
                    }
                    cout << endl;
                    throw runtime_error("Duplicate gene name in hdf5 file " + fileName);
                }
            }
        }

        // Read the cell names.
        const H5::DataSet cellNamesDataSet = file.openDataSet(groupName + "barcodes");
        vector<string> hdf5CellNames;
        hdf5::read(cellNamesDataSet, hdf5CellNames);
        cout << "Found " << hdf5CellNames.size() << " cells." << endl;

        // Read the index pointers.
        // These can be used to locate the information for each cell in the
        // data and indices vectors.
        const H5::DataSet indexPointersDataSet = file.openDataSet(groupName + "indptr");
        vector<uint64_t> indexPointers;
        hdf5::read(indexPointersDataSet, indexPointers);
        if (indexPointers.size() != hdf5CellNames.size() + 1) {
            const string message =
                "Unexpected length of index pointers in hdf5 file " + fileName +
                ": " + lexical_cast<string>(indexPointers.size() +
                ". Expected " + lexical_cast<string>(hdf5CellNames.size() + 1) + ".");
        }



        // Add the genes.
        // We want to add them independently of the cells, so they all get added, even the ones
        // for which all cells have zero count.
        for (const string& hdf5GeneName : hdf5GeneNames) {
            addGene(hdf5GeneName);
        }



        // Main loop over the cells.
        // For each cell we read the appropriate portion of the data and indices vector.
        // Data contains expression counts and indices contain gene indices into the hdf5GeneNames vector.
        // Note that the gene indices do not necessarily agree with the GeneId as stored in the
        // Expression Matrix (this is in particular true if some genes were already defined
        // before the call to addCellsFromHdf5).
        // Similarly, the index i of this loop is not necessarily the same as the CellId
        // of the cell being added.
        // The data and indices vectors will be used to hold the data and indices for
        // a single cell.
        // The metaData and expressionCounts vector are used as arguments to addCells.
        // The four vectors are defined here to avoid reallocation inside the loop.
        vector<uint32_t> data;
        vector<uint64_t> indices;	// Silly, but that's the way 10X does it.
        vector<pair<string, string> > metaData(1, make_pair("CellName", ""));
        vector<pair<string, float> > expressionCounts;
        const H5::DataSet dataDataSet = file.openDataSet(groupName + "data");
        const H5::DataSet indicesDataSet = file.openDataSet(groupName + "indices");
        for (size_t i = 0; i < hdf5CellNames.size(); i++) {
            if ((i % 1000) == 0) {
                cout << timestamp << "Added " << i << " cells of " << hdf5CellNames.size() << endl;
            }

            // Read the data and indices for this cell.
            const uint64_t offset = indexPointers[i];
            const uint64_t n = indexPointers[i + 1] - offset;
            hdf5::read(dataDataSet, offset, n, data);
            hdf5::read(indicesDataSet, offset, n, indices);
            CZI_ASSERT(data.size() == n);
            CZI_ASSERT(indices.size() == n);

            // cout << "offset = " << offset << endl;
            // cout << "n = " << n << endl;

            // Add this cell.
            expressionCounts.clear();
            metaData.front().second = hdf5CellNames[i];
            for (size_t j = 0; j < n; j++) {
                // cout << j << " " << indices[j] << endl;
                const uint64_t hdf5GeneIndex = indices[j];
                CZI_ASSERT(hdf5GeneIndex < hdf5GeneNames.size());
                expressionCounts.push_back(make_pair(hdf5GeneNames[hdf5GeneIndex], float(data[j])));
            }
            try {
                addCell(metaData, expressionCounts);
            } catch (...) {
                cout << "Error occurred adding cell " << i << " ";
                cout << hdf5CellNames[i] << " from HDF5 file " << fileName << endl;
                cout << "Index pointers information: " << indexPointers[i];
                cout << " " << indexPointers[i + 1] << " " << n << endl;
                cout << ". Details follow." << endl;
                for (size_t j = 0; j < n; j++) {
                    cout << j << " " << indices[j] << " " << hdf5GeneNames[indices[j]] << " " << data[j] << endl;
                }
                cout << "Error occurred adding cell " << i << " ";
                cout << hdf5CellNames[i] << " from HDF5 file " << fileName  << endl;
                cout << "Index pointers information: " << indexPointers[i] << " ";
                cout << indexPointers[i + 1] << " " << n << endl;
                cout << "See details above." << endl;
                throw;
            }
        }



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

    cout << "There are " << cellCount() << " cells and " << geneCount() << " genes." << endl;
}
