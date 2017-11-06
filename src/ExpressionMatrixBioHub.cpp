/*******************************************************************************

Functionality to read expression matrix data from the BioHub pipeline.



JULY 2017, ILLUMINA DATA: addCellsFromBioHub1

This version of the BioHub pipeline creates three files for each plate:

- A csv file containing expression counts by cell, with one row per cell
  and one column per gene plus a header line containing gene names.
  The name of this file is specified as the first parameter to addCellsFromBioHub.
  The first column contains the cell name.
  Additional columns are also allowed before and after the last gene.
  These columns are treated as per-cell meta data.
  The number of these initial and final columns is specified in the second and
  third parameter to addCellsFromBioHub.
- A csv file containing plate meta data. The first row contains plate names.
  The plate name is obtained from the name of the expression counts file
  (it is the portion in the file name that preceded the first period).
  The row in this file corresponding to that plate name is used to assign meta
  data to all the cells.
- A csv file containing cell meta data, with one row for each cell and one
  column for each meta data field. The first column contains the cell name
  (which matches the cell name used in the expression count file).
  The cells in this file are not required
  to be in the same order as the cells in the expression counts file.

The names of the first two files are passed as arguments to addCellGFromBioHub.
The additional cell meta data in the third file can then be added using a call
to addCellMetaData.



SEPTEMBER 2017, 10X GENOMICS DATA: addCellsFromBioHub2

This uses a master file which is a comma separated file containing
a header line containing column names plus one line for each plate.
The name of this master file is given as the only argument
to addCellsFromBioHub2.
A single call to addCellsFromBioHub2 is used to load all of the cells.

The following two columns are mandatory:
    - PlateName: contains the name of the plate that
      each line in the file refers to.
    - aws-s3-hdf5: contains the AWS S3 path name of the hdf5 for the plate.
      A valid path must begin with "s3://" followed by the name of the
      S3 bucket that contains the data.

The remaining columns are treated as plate meta data.
The meta data for a plate are added to all of the cells found on the plate.

Cells are named using the pattern PlateName-Barcode, where the barcode
is as obtained from the HDF5 file for each plate.



NOVEMBER 2017, ILLUMINA DATA: addCellsFromBioHub3

For each of the plates to be processed, there is an expression counts file
and a cell meta data file. In both files, each column corresponds to a cell,
in arbitrary order.
Only cells that are present in both files are added.
Both files are tab separated, and the header line has
one fewer fields than all remaining lines.

Arguments:

plateListFileName:
Name of the file containing the list of plate names to be processed, one per line.

plateMetaDataFileName
Name of the csv file containing per-plate meta data.
The first column of each line contains a plate name.
This meta data for each plate is added to all the cells of the plate.
This file can contain data for a superset of the plate names
specified in plateListFileName. Only the plates
specified in plateListFileName will be processed.

expressionCountsFileNameSuffix
Suffix to be added to a plate name to obtain the name of the corresponding expression counts file.

metaDataFileNameSuffix
Suffix to be added to a plate name to obtain the name of the corresponding per-cell meta data file.


*******************************************************************************/



#include "ExpressionMatrix.hpp"
#include "Aws.hpp"
#include "timestamp.hpp"
#include "tokenize.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "fstream.hpp"
#include "iterator.hpp"



//JULY 2017, ILLUMINA DATA: addCellsFromBioHub1
// See the beginning of this file for more information.
void ExpressionMatrix::addCellsFromBioHub1(
    const string& expressionCountsFileName, // The name of the csv file containing expression counts.
    size_t initialMetaDataCount,            // The number of initial columns containing meta data.
    size_t finalMetaDataCount,              // The number of final columns containing meta data.
    const string& plateMetaDataFileName     // The name of the file containing per-plate meta data.
    )
{

    // Open the expression counts file.
    ifstream expressionCountsFile(expressionCountsFileName);
    if(!expressionCountsFile) {
        throw runtime_error("Error opening " + expressionCountsFileName);
    }

    // Read the first row.
    string line;
    getline(expressionCountsFile, line);
    if(!expressionCountsFile) {
        throw runtime_error("Error reading the header line from file " + expressionCountsFileName);
    }

    // Parse the first row.
    vector<string> tokens;
    tokenize(",", line, tokens, true);
    const size_t minimumExpectedTokenCount = initialMetaDataCount + finalMetaDataCount + 2;
    if(tokens.size() < minimumExpectedTokenCount) {
        throw runtime_error("Insufficient number of tokens in first line of file " + expressionCountsFileName +
            " Expected at least " + lexical_cast<string>(minimumExpectedTokenCount) +
            "tokens, found " + lexical_cast<string>(tokens.size()) + "."
            );
    }
    const size_t tokenCount = tokens.size();
    const size_t genesInCsvFileCount = tokenCount - 1 - initialMetaDataCount - finalMetaDataCount;

    // Compute the indexes of the begin/end token for initial meta data, expression counts, final meta data.
    const size_t initialMetaDataBegin = 1;
    const size_t initialMetaDataEnd = initialMetaDataBegin + initialMetaDataCount;
    const size_t expressionCountBegin = initialMetaDataEnd;
    const size_t expressionCountEnd = expressionCountBegin + genesInCsvFileCount;
    const size_t finalMetaDataBegin = expressionCountEnd;
    const size_t finalMetaDataEnd = finalMetaDataBegin + finalMetaDataCount;
    CZI_ASSERT(finalMetaDataEnd == tokens.size());

    // Initial meta data names.
    const vector<string> initialMetaDataNames(tokens.begin() + initialMetaDataBegin,
        tokens.begin() + initialMetaDataEnd);
    CZI_ASSERT(initialMetaDataNames.size() == initialMetaDataCount);

    // Names of the genes in this csv file.
    const vector<string> geneNamesInCsvFile(tokens.begin() + expressionCountBegin, tokens.begin() + expressionCountEnd);
    CZI_ASSERT(geneNamesInCsvFile.size() == genesInCsvFileCount);

    // Initial meta data names.
    const vector<string> finalMetaDataNames(tokens.begin() + finalMetaDataBegin, tokens.begin() + finalMetaDataEnd);
    CZI_ASSERT(finalMetaDataNames.size() == finalMetaDataCount);

    // Some messages.
    cout << timestamp << "Working on file " << expressionCountsFileName << endl;
    cout << "This file  contains expression counts for ";
    cout << genesInCsvFileCount << " genes plus" << endl;
    cout << "cell names and " << initialMetaDataCount << " initial columns and ";
    cout << finalMetaDataCount << " final columns of cell meta data" << endl;
    cout << "for a total " << tokens.size() << " tokens per line." << endl;
    cout << "First gene in this csv file: " << geneNamesInCsvFile.front() << endl;
    cout << "Last gene in this csv file: " << geneNamesInCsvFile.back() << endl;
    /*
     cout << "Initial meta data names:\n" << endl;
     copy(initialMetaDataNames.begin(), initialMetaDataNames.end(), ostream_iterator<string>(cout, "\n"));
     cout << "Final meta data names: ";
     copy(finalMetaDataNames.begin(), finalMetaDataNames.end(), ostream_iterator<string>(cout, "\n"));
     cout << flush;
     */

    // Add all the genes. We add them here so they all get added, even those that have
    // zero counts for all cells.
    for(const string& geneName : geneNamesInCsvFile) {
        addGene(geneName);
    }

    // Extract the plate name from the name of the expression counts file.
    // It is the portion before the period in the file name.
    const string expressionCountsFileNameOnly = boost::filesystem::path(expressionCountsFileName).filename().string();
    const string plateName = expressionCountsFileNameOnly.substr(0, expressionCountsFileNameOnly.find_first_of('.'));
    cout << "Plate name is " << plateName << endl;

    // Get the meta data for the plate corresponding to this expression counts file.
    // This meta data will be assigned to all of the cells in the file.
    vector<pair<string, string> > plateMetaData;
    getPlateMetaDataFromBioHub(plateName, plateMetaDataFileName, plateMetaData);



    // Vectors to contain expression counts and meta data for a single cell.
    vector<pair<string, float> > cellExpressionCounts;
    vector<pair<string, string> > cellMetaData;
    cellMetaData.push_back(make_pair("CellName", ""));
    cellMetaData.insert(cellMetaData.end(), plateMetaData.begin(), plateMetaData.end());
    for(const string& metaDataName : initialMetaDataNames) {
        cellMetaData.push_back(make_pair(metaDataName, ""));
    }
    for(const string& metaDataName : finalMetaDataNames) {
        cellMetaData.push_back(make_pair(metaDataName, ""));
    }



	// Read the cells, one per line.
    size_t newCellCount = 0;
    while(true) {

        // Get a line and parse it.
        getline(expressionCountsFile, line);
        if(!expressionCountsFile) {
            break;
        }
        tokenize(",", line, tokens, true);
        if(tokens.size() != tokenCount) {
            cout << line << endl;
            cout << "The previous line has " << tokens.size() << " tokens. Expected " << tokenCount << endl;
            throw runtime_error("Invalid number of tokens in line of input expression count file.");
        }

        // Store the initial and final meta data names.
        cellMetaData.front().second = tokens.front();
        for(size_t i = 0; i != initialMetaDataCount; i++) {
            cellMetaData[i + 1 + plateMetaData.size()].second = tokens[initialMetaDataBegin + i];
        }
        for(size_t i = 0; i != finalMetaDataCount; i++) {
            cellMetaData[i + 1 + plateMetaData.size() + initialMetaDataCount].second = tokens[finalMetaDataBegin + i];
        }

        // Store the expression counts.
        cellExpressionCounts.clear();
        for(size_t i = 0; i < genesInCsvFileCount; i++) {
            const string& token = tokens[expressionCountBegin + i];
            float count;
            try {
                count = lexical_cast<float>(token);
            } catch (bad_lexical_cast) {
                throw runtime_error("Invalid format of expression count " + token + " for cell " + tokens.front());
            }
            if(count == 0.) {
                continue;
            }
            if(count > 0.) {
                cellExpressionCounts.push_back(make_pair(geneNamesInCsvFile[i], count));
            } else {
                throw runtime_error("Negative expression count " + token + " for cell " + tokens.front());
            }
        }

        // Add the cell.
        addCell(cellMetaData, cellExpressionCounts);
        ++newCellCount;
    }
    cout << "Read expression counts for " << newCellCount << " cells." << endl;
    cout << "There are " << cellCount() << " cells and " << geneCount() << " genes." << endl;


}



void ExpressionMatrix::getPlateMetaDataFromBioHub(
    const string& plateName,
    const string&plateMetaDataFileName,
    vector<pair<string, string> >& plateMetaData)
{
    // Open the plate meta data file.
    ifstream plateMetaDataFile(plateMetaDataFileName);
    if(!plateMetaDataFile) {
        throw runtime_error("Error opening " + plateMetaDataFileName);
    }

    // Read the first row.
    string line;
    getline(plateMetaDataFile, line);
    if(!plateMetaDataFile) {
        throw runtime_error("Error reading the header line from file " + plateMetaDataFileName);
    }

    // Parse the first row. They are the meta data names.
    vector<string> metaDataNames;
    tokenize(",", line, metaDataNames, true);
    if(metaDataNames.size() == 0) {
        throw runtime_error("Invalid format of plate meta data file " + plateMetaDataFileName);
    }
    metaDataNames.front() = "PlateName";



	// Find the row for the plate we are looking for.
    vector<string> metaDataValues;
    while(true) {

        getline(plateMetaDataFile, line);
        if(!plateMetaDataFile) {
            break;
        }
        tokenize(",", line, metaDataValues, true);

        if(metaDataValues.size() != metaDataNames.size()) {
            throw runtime_error("Unexpected number of meta data values for plate " + plateName);
        }

        if(metaDataValues.front() != plateName) {
            continue;
        }

        // We found the line for this plate.
        // Now we have both names and values of the meta data for this cell.
        plateMetaData.clear();
        for(size_t i = 0; i < metaDataNames.size(); i++) {
            plateMetaData.push_back(make_pair(metaDataNames[i], metaDataValues[i]));
        }
        return;

    }

    throw runtime_error(
        "Did not find line for plate " + plateName +
            " in plate meta data file " + plateMetaDataFileName
            );
}



void ExpressionMatrix::addCellMetaData(const string& cellMetaDataFileName)
{
    // Open the cell meta data file.
    ifstream cellMetaDataFile(cellMetaDataFileName);
    if(!cellMetaDataFile) {
        throw runtime_error("Error opening " + cellMetaDataFileName);
    }
    cout << timestamp << "Adding cell meta data in " << cellMetaDataFileName << endl;

    // Read the first row.
    string line;
    getline(cellMetaDataFile, line);
    if(!cellMetaDataFile) {
        throw runtime_error("Error reading the header line from file " + cellMetaDataFileName);
    }

    // Parse the first row. They are the meta data names (except for the first one which is ignored).
    vector<string> metaDataNames;
    tokenize(",", line, metaDataNames, true);
    if(metaDataNames.size() == 0) {
        throw runtime_error("Invalid format of cell meta data file " + cellMetaDataFileName);
    }



    // Process each row.
    vector<string> metaDataValues;
    while(true) {

        getline(cellMetaDataFile, line);
        if(!cellMetaDataFile) {
            break;
        }
        tokenize(",", line, metaDataValues, true);

        if(metaDataValues.size() != metaDataNames.size()) {
            throw runtime_error("Unexpected number of meta data values in " + cellMetaDataFileName);
        }

        // Locate the cell.
        const string& cellName = metaDataValues.front();
        const StringId cellNameStringId = cellNames(cellName);
        if(cellNameStringId == cellNames.invalidStringId) {
            throw runtime_error("Attempt to add meta data for undefined cell " + cellName);
        }
        const CellId cellId = CellId(cellNameStringId);
        CZI_ASSERT(cellId < cellCount());

        for(size_t i = 1; i < metaDataNames.size(); i++) {
            if(metaDataNames[i].size() > 0) {
                setCellMetaData(cellId, metaDataNames[i], metaDataValues[i]);
            }
        }
    }

}



// SEPTEMBER 2017, 10X GENOMICS DATA: addCellsFromBioHub2
// See the beginning of this file for more information.
void ExpressionMatrix::addCellsFromBioHub2(
    const string& plateFileName,
    double totalExpressionCountThreshold
    )
{
    // Open the plate file.
    ifstream plateFile(plateFileName);
    if(!plateFile) {
        throw runtime_error("Plate file " + plateFileName + " not found or could not be opened.");
    }

    // Read the header line.
    string line;
    getline(plateFile, line);
    if(!plateFile || line.size()==0) {
        throw runtime_error("Unable to read header line from plate file " + plateFileName);
    }

    // Parse the header line to get the column names.
    if(line[line.size()-1] == 13) {
        line.resize(line.size()-1); // Remove Windows style line end if necessary.
    }
    vector<string> columnNames;
    tokenize(",", line, columnNames);

    // Find the index of the column containing the PlateName.
    const auto plateNameIt = find(columnNames.begin(), columnNames.end(), "PlateName");
    if(plateNameIt == columnNames.end()) {
        throw runtime_error("Plate file " + plateFileName + " does not have PlateName column.");
    }
    const size_t plateNamePosition = plateNameIt - columnNames.begin();

    // Find the index of the column containing the aws-s3-hdf5 path.
    const auto awsPathIt = find(columnNames.begin(), columnNames.end(), "aws-s3-hdf5");
    if(awsPathIt == columnNames.end()) {
        throw runtime_error("Plate file " + plateFileName + " does not have aws-s3-hdf5 column.");
    }
    const size_t awsPathPosition = awsPathIt - columnNames.begin();


    // Prepare the vector of cell meta data.
    vector< pair<string, string> > cellMetaData;
    for(const string& columnName: columnNames) {
        cellMetaData.push_back(make_pair(columnName, ""));
    }



    // Main loop over the remaining lines of the plate file.
    // Each line corresponds to a plate.
    size_t plateCount = 0;
    for(; ; ++plateCount) {

        // Read this line.
        getline(plateFile, line);
        if(!plateFile | line.empty()) {
            break;
        }
        if(line[line.size()-1] == 13) {
            line.resize(line.size()-1); // Remove Windows style line end if necessary.
        }

        // Parse it.
        vector<string> columnValues;
        tokenize(",", line, columnValues);
        if(columnValues.size() != columnNames.size()) {
            cout << line << endl;
            cout << "Found " << columnValues.size() << " tokens in the above line, expected ";
            cout << columnValues.size() << endl;
            throw runtime_error("Unexpected number of tokens in line of plate file " + plateFileName);

        }

        // Fill in the cell meta data for this plate.
        for(size_t i=0; i<columnNames.size(); i++) {
            cellMetaData[i].second = columnValues[i];
        }

        // Extract the PlateName and the aws path.
        CZI_ASSERT(plateNamePosition < columnValues.size());
        CZI_ASSERT(awsPathPosition < columnValues.size());
        const string& plateName = columnValues[plateNamePosition];
        const string awsPath = columnValues[awsPathPosition];
        cout << timestamp << " Working on plate " << plateCount << " " << plateName << endl;

        // Make a local copy of the hdf5 file.
        const string uuid = boost::uuids::to_string(boost::uuids::uuid(boost::uuids::random_generator()()));
        const string localPath = string("/dev/shm/aws-") + uuid;
        const string command = "aws s3 cp --quiet " + awsPath + " " + localPath;
        const int returnCode = ::system(command.c_str());
        if(returnCode!=0) {
            throw runtime_error("Error " +
                lexical_cast<string>(returnCode) +
                " from command: " + command);
        }

        // Add cells from this hdf5 file.
        addCellsFromHdf5(localPath, plateName, cellMetaData, totalExpressionCountThreshold);


        // Remove the local copy.
        boost::filesystem::remove(localPath);

    }
    cout << timestamp << "Processed " << plateCount << " plates." << endl;
}



// November 2017, Illumina data.
// See the beginning of this file for more information.
void ExpressionMatrix::addCellsFromBioHub3(
    const string& expressionCountsFileName,             // The name of the csv file containing expression counts.
    const string& expressionCountsFileSeparators,
    const vector<pair<string, string> >& plateMetaDataWithoutCellName  // Meta data that will be added to all cells.
    )
{
    // Create cell meta data including a slot for the cell name.
    vector<pair<string, string> > plateMetaData;
    plateMetaData.push_back(make_pair("CellName", ""));
    copy(plateMetaDataWithoutCellName.begin(), plateMetaDataWithoutCellName.end(),
        back_inserter(plateMetaData));

    // See how many cells we have in this plate.
    const size_t cellCountInPlate = countTokensInSecondLine(expressionCountsFileName, expressionCountsFileSeparators) - 1;


    // Open the expression counts file for this plate.
    ifstream expressionCountsFile = ifstream(expressionCountsFileName);
    if(!expressionCountsFile) {
        throw runtime_error("Error opening " + expressionCountsFileName);
    }


    // Read the cell names from the first line.
    string line;
    getline(expressionCountsFile, line);
    removeWindowsLineEnd(line);
    vector<string> cellNamesInPlate;
    tokenize(expressionCountsFileSeparators, line, cellNamesInPlate);
    if(cellNamesInPlate.size() == cellCountInPlate) {
        // Do nothing.
    } else if(cellNamesInPlate.size() == cellCountInPlate+1) {
        cellNamesInPlate.erase(cellNamesInPlate.begin());
    } else {
        throw runtime_error(
            "File " + expressionCountsFileName +
            " has inconsistent number of tokens in first two lines.");
    }



    // Read the expression file to create expression vectors for all the cells in this plate.
    vector< vector< pair<string, float> > > cellsInPlateExpressionVectors(cellCountInPlate);
    vector<string> tokens;
    for(size_t lineNumber=2; ; ++lineNumber) {

        // Get a line.
        getline(expressionCountsFile, line);
        if(!expressionCountsFile) {
            break;
        }
        removeWindowsLineEnd(line);
        tokenize(expressionCountsFileSeparators, line, tokens);
        if(tokens.size() != cellCountInPlate+1) {
            throw runtime_error(
                "Invalid number of tokens in file "
                + expressionCountsFileName +
                " line " + lexical_cast<string>(lineNumber));
        }

        // Add the gene name.
        const string& geneName = tokens.front();

        // Add all the non-zero expression counts.
        for(size_t i=0; i<cellCountInPlate; i++) {
            const string& token = tokens[i+1];
            if(token=="0" || token=="0.") {
                continue;
            }
            float expressionCount;
            try {
                expressionCount = lexical_cast<float>(token);
            } catch(bad_lexical_cast) {
                throw runtime_error("Invalid expression count at line " +
                    lexical_cast<string>(lineNumber) + " of file " + expressionCountsFileName);
            }
            cellsInPlateExpressionVectors[i].push_back(make_pair(geneName, expressionCount));
        }

    }


    // Now we can add the cells.
    for(size_t i=0; i<cellCountInPlate; i++) {
        plateMetaData.front().second = cellNamesInPlate[i];
        addCell(plateMetaData, cellsInPlateExpressionVectors[i]);
    }




}
