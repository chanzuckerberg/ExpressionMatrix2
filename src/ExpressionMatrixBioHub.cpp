/*******************************************************************************

Functionality to read expression matrix data from the BioHub pipeline.

The BioHub pipeline creates three files for each plate:

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
  column for each meta data field. The first column containes the cell name
  (which matches the cell name used in the excpression count file).
  The cells in this file are not required
  to be in the same order as the cells in the expression counts file.

The names of the first two files are passed as arguments to addCellGFromBioHub.
The additional cell meta datga in the third file can then be added using a call
to addCellMetaData.

*******************************************************************************/



#include "ExpressionMatrix.hpp"
#include "timestamp.hpp"
#include "tokenize.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/filesystem.hpp>

#include "fstream.hpp"
#include "iterator.hpp"



void ExpressionMatrix::addCellsFromBioHub(
	const string& expressionCountsFileName,	// The name of the csv file containing expression counts.
	size_t initialMetaDataCount,			// The number of initial columns containing meta data.
	size_t finalMetaDataCount,				// The number of final columns containing meta data.
	const string& plateMetaDataFileName 	// The name of the file containing per-plate meta data.
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
	const vector<string> initialMetaDataNames(tokens.begin()+initialMetaDataBegin, tokens.begin()+initialMetaDataEnd);
	CZI_ASSERT(initialMetaDataNames.size() == initialMetaDataCount);

	// Names of the genes in this csv file.
	const vector<string> geneNamesInCsvFile(tokens.begin()+expressionCountBegin, tokens.begin()+expressionCountEnd);
	CZI_ASSERT(geneNamesInCsvFile.size() == genesInCsvFileCount);

	// Initial meta data names.
	const vector<string> finalMetaDataNames(tokens.begin()+finalMetaDataBegin, tokens.begin()+finalMetaDataEnd);
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
	for(const string& geneName: geneNamesInCsvFile) {
		addGene(geneName);
	}


	// Extract the plate name from the name of the expression counts file.
	// It is the portion before the period in the file name.
	const string expressionCountsFileNameOnly = boost::filesystem::path(expressionCountsFileName).filename().string();
	const string plateName = expressionCountsFileNameOnly.substr(0, expressionCountsFileNameOnly.find_first_of('.'));
	cout << "Plate name is " << plateName << endl;

	// Get the meta data for the plate corresponding to this expression counts file.
	// This meta data will be assigned to all of the cells in the file.
	vector< pair<string, string> > plateMetaData;
	getPlateMetaDataFromBioHub(plateName, plateMetaDataFileName, plateMetaData);



	// Vectors to contain expression counts and meta data for a single cell.
    vector< pair<string, float> > cellExpressionCounts;
    vector< pair<string, string> > cellMetaData;
    cellMetaData.push_back(make_pair("CellName", ""));
    cellMetaData.insert(cellMetaData.end(), plateMetaData.begin(), plateMetaData.end());
    for(const string& metaDataName: initialMetaDataNames) {
    	cellMetaData.push_back(make_pair(metaDataName, ""));
    }
    for(const string& metaDataName: finalMetaDataNames) {
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
		for(size_t i=0; i!=initialMetaDataCount; i++) {
			cellMetaData[i + 1 + plateMetaData.size()].second = tokens[initialMetaDataBegin+i];
		}
		for(size_t i=0; i!=finalMetaDataCount; i++) {
			cellMetaData[i + 1 + plateMetaData.size() + initialMetaDataCount].second = tokens[finalMetaDataBegin+i];
		}

		// Store the expression counts.
		cellExpressionCounts.clear();
		for(size_t i=0; i<genesInCsvFileCount; i++) {
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
	vector< pair<string, string> >& plateMetaData)
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
		for(size_t i=0; i<metaDataNames.size(); i++) {
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
