#include "ExpressionMatrix.hpp"
#include "CellSimilarityGraph.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
#include "tokenize.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/regex.hpp>

#include "fstream.hpp"
#include "iostream.hpp"
#include "utility.hpp"
#include "vector.hpp"
#include <numeric>
#include <sstream>



// Construct a new expression matrix. All binary data for the new expression matrix
// will be stored in the specified directory. If the directory does not exist,
// it will be created. If the directory already exists, any previous
// expression matrix stored in the directory will be overwritten by the new one.
ExpressionMatrix::ExpressionMatrix(
    const string& directoryName,
    const ExpressionMatrixCreationParameters& parameters) :
    directoryName(directoryName)
{
    // If directory already exists, don't do anything.
	// This ensures that we don't delete or overwrite anything,
	// and that the directory does not contain stale data.
    if(boost::filesystem::exists(directoryName)) {
    	throw runtime_error("Directory "  + directoryName + " already exists.");
    }

    // Create the directory. This guarantees that we start with an empty directory.
    if(!boost::filesystem::create_directory(directoryName)) {
    	throw runtime_error("Unable to create directory "  + directoryName);
    }

    geneNames.createNew(directoryName + "/" + "GeneNames", parameters.geneCapacity);
    cells.createNew(directoryName + "/" + "Cells");
    cellNames.createNew(directoryName + "/" + "CellNames", parameters.cellCapacity);
    cellMetaData.createNew(directoryName + "/" + "CellMetaData");
    cellMetaDataNames.createNew(directoryName + "/" + "CellMetaDataNames", parameters.cellMetaDataNameCapacity);
    cellMetaDataValues.createNew(directoryName + "/" + "CellMetaDataValues", parameters.cellMetaDataValueCapacity);
    cellMetaDataNamesUsageCount.createNew(directoryName + "/" + "CellMetaDataNamesUsageCount");
    cellExpressionCounts.createNew(directoryName + "/" + "CellExpressionCounts");
    largeCellExpressionCounts.createNew(directoryName + "/" + "LargeCellExpressionCounts");

    // Initialize the CellSets.
    cellSets.createNew(directoryName);
    vector<CellId> emptyCellSet;
    cellSets.addCellSet("AllCells", emptyCellSet);

    // Initialize the gene sets.
    geneSets["AllGenes"].createNew(directoryName + "/GeneSet-AllGenes");

    // Sanity checks.
    CZI_ASSERT(cellNames.size() == cells.size());
    CZI_ASSERT(cellMetaData.size() == cells.size());
    CZI_ASSERT(cellExpressionCounts.size() == cells.size());
    CZI_ASSERT(cellSets.cellSets["AllCells"]->size() == cells.size());

    // Fill the table containing commands known to the http server.
    fillServerFunctionTable();
}



// Access a previously created expression matrix stored in the specified directory.
ExpressionMatrix::ExpressionMatrix(const string& directoryName) :
    directoryName(directoryName)
{
    // Access the binary data with read-write access, so we can add new cells
    // and perform other operations that change the state on disk.
    geneNames.accessExistingReadWrite(directoryName + "/" + "GeneNames");
    cells.accessExistingReadWrite(directoryName + "/" + "Cells");
    cellNames.accessExistingReadWrite(directoryName + "/" + "CellNames");
    cellMetaData.accessExistingReadWrite(directoryName + "/" + "CellMetaData");
    cellMetaDataNames.accessExistingReadWrite(directoryName + "/" + "CellMetaDataNames");
    cellMetaDataValues.accessExistingReadWrite(directoryName + "/" + "CellMetaDataValues");
    cellMetaDataNamesUsageCount.accessExistingReadWrite(directoryName + "/" + "CellMetaDataNamesUsageCount");
    cellExpressionCounts.accessExistingReadWrite(directoryName + "/" + "CellExpressionCounts");
    largeCellExpressionCounts.accessExistingReadWrite(directoryName + "/" + "LargeCellExpressionCounts");
    cellSets.accessExisting(directoryName);



    // Access the gene sets.
	using boost::filesystem::directory_iterator;
	boost::regex regex(directoryName + "/GeneSet-(.*)-Genes");
	for(auto it=directory_iterator(directoryName); it!=directory_iterator(); ++it) {
		const string fileName = it->path().string();
		boost::smatch regexMatchResults;
		if(!boost::regex_match(fileName, regexMatchResults, regex)) {
			continue;
		}
		CZI_ASSERT(regexMatchResults.size() == 2);
		const string& geneSetName = regexMatchResults[1];
		geneSets[geneSetName].accessExisting(directoryName + "/GeneSet-" + geneSetName);
	}



    // Sanity checks.
    CZI_ASSERT(cellNames.size() == cells.size());
    CZI_ASSERT(cellMetaData.size() == cells.size());
    CZI_ASSERT(cellExpressionCounts.size() == cells.size());
    CZI_ASSERT(cellSets.cellSets["AllCells"]->size() == cells.size());
    CZI_ASSERT(cellMetaDataNamesUsageCount.size() == cellMetaDataNames.size());
    CZI_ASSERT(geneSets["AllGenes"].size() == geneCount());

    // Fill the table containing commands known to the http server.
    fillServerFunctionTable();
}



// Add a gene.
// Returns true if the gene was added, false if it was already present.
bool ExpressionMatrix::addGene(const string& geneName)
{
    const StringId stringId = geneNames(geneName);
    if(stringId == geneNames.invalidStringId) {
		const GeneId geneId = GeneId(geneNames[geneName]);
		geneSets["AllGenes"].addGene(geneId);
		return true;
    } else {
    	return false;	// Was already present.
    }
}



// Add a cell to the expression matrix.
// The meta data is passed as a vector of names and values, which are all strings.
// The cell name should be entered as meta data "CellName".
// The expression counts for each gene are passed as a vector of pairs
// (gene names, count).
// Returns the id assigned to this cell.
// This changes the metaData vector so the CellName entry is the first entry.
// It also changes the expression counts - it sorts them by decreasing count.
CellId ExpressionMatrix::addCell(
    vector< pair<string, string> >& metaData,
    vector< pair<string, float> >& expressionCounts,
    size_t maxTermCountForApproximateSimilarityComputation)
{
    // Check that we don't overflow the CellId type.
    CZI_ASSERT(CellId(cells.size()) < std::numeric_limits<CellId>::max());

    // Make sure the cellName entry exists and place it at the beginning
    // of the meta data.
    bool cellNameWasFound = false;
    for(auto& p: metaData) {
        if(p.first == "CellName") {
            cellNameWasFound = true;
            swap(p, metaData.front());
            break;
        }
    }
    if(!cellNameWasFound) {
    	cout << "CellName is missing from the following cell meta data:" << endl;
    	for(const auto& p: metaData) {
    		cout << p.first << " " << p.second << endl;
    	}
        throw std::runtime_error("CellName missing from cell meta data.");
    }
    CZI_ASSERT(metaData.front().first == "CellName");

    // Sort the rest of the meta data.
    // sort(metaData.begin()+1, metaData.end());

    // Check that we don't already have this cell name.
    const string& cellName = metaData.front().second;
    StringId cellNameStringId = cellNames(cellName);
    if(cellNameStringId != invalidCellId) {
        throw runtime_error("Cell name " + cellName + " already exists.");
    }

    // Store the cell name.
    cellNameStringId = cellNames[cellName];



    // Store the cell meta data.
    cellMetaData.push_back();
    for(const auto& p: metaData) {

        // Get the StringId for the name.
        const StringId nameId = cellMetaDataNames[p.first];

        // Increment the usage count for this name.
        incrementCellMetaDataNameUsageCount(nameId);

        // Get the StringId for the value.
        const StringId valueId = cellMetaDataValues[p.second];

        // Store the (name,value) pair.
        cellMetaData.push_back(make_pair(nameId, valueId));
    }


    // Store the expression counts.
    Cell cell;
    cell.sum1 = 0.;
    cell.sum2 = 0.;
    cellExpressionCounts.appendVector();
    for(const auto& p: expressionCounts) {
        const StringId geneId = geneNames[p.first];
        const float value = p.second;
        if(value < 0) {
            throw runtime_error("Negative expression count encountered.");
        }
        cell.sum1 += value;
        cell.sum2 += value*value;
        cellExpressionCounts.append(make_pair(geneId, value));
    }
    cell.norm2 = sqrt(cell.sum2);
    cell.norm1Inverse = 1./cell.norm1();
    cell.norm2Inverse = 1./cell.norm2;

    // Sort the expression counts we just stored by GeneId.
    const auto storedExpressionCounts = cellExpressionCounts[cellExpressionCounts.size()-1];
    sort(storedExpressionCounts.begin(), storedExpressionCounts.end());



    // Verify that all the gene ids in the expression counts we just stored are distinct.
    for(size_t i=1; i<storedExpressionCounts.size(); i++) {
    	if(storedExpressionCounts[i-1].first == storedExpressionCounts[i].first) {
    		const string duplicateGeneName = geneNames[storedExpressionCounts[i].first];
    		throw runtime_error("Duplicate expression count for cell " + cellName + " gene " + duplicateGeneName);
    	}
    }

    // We need to sort the input expression counts by decreasing count.
    sort(expressionCounts.begin(), expressionCounts.end(),
        OrderPairsBySecondGreaterThenByFirstLess< pair<string, float> >());

    // Store the maxTermCountForApproximateSimilarityComputation largest expression counts
    // for use by computeApproximateCellSimilarity.
    const size_t numberToKeep = min(expressionCounts.size(), maxTermCountForApproximateSimilarityComputation);
    largeCellExpressionCounts.appendVector();
    for(size_t i=0; i<numberToKeep; i++) {
        const pair<string, float>& p = expressionCounts[i];
        const StringId geneId = geneNames[p.first];
        const float value = p.second;
        largeCellExpressionCounts.append(make_pair(geneId, value));
    }

    // Sort the large expression counts we just stored by GeneId.
    const auto storedLargeExpressionCounts = largeCellExpressionCounts[cellExpressionCounts.size()-1];
    sort(storedLargeExpressionCounts.begin(), storedLargeExpressionCounts.end());

    // Store cell.sum1LargeExpressionCounts and cell.sum2LargeExpressionCounts
    // for use by computeApproximateCellSimilarity.
    cell.sum1LargeExpressionCounts = 0.;
    cell.sum2LargeExpressionCounts = 0.;
    for(const auto& p: storedLargeExpressionCounts) {
        const double& count = p.second;
        cell.sum1LargeExpressionCounts += count;
        cell.sum2LargeExpressionCounts += count*count;
    }


#if 0

    // Additional code used to compute bounds on approximate cell similarity.
    // None of the bounds turned out to be sufficiently tight, so we are not using these.

    cell.largeExpressionCountsNorm2 = sqrt(cell.sum2LargeExpressionCounts);
    double sum2SmallExpressionCount = 0.;
    for(size_t i=numberToKeep; i<expressionCounts.size(); i++) {
        const auto& p = expressionCounts[i];
        const double& count = p.second;
        sum2SmallExpressionCount += count*count;
    }
    cell.smallExpressionCountsNorm2 = sqrt(sum2SmallExpressionCount);

    // L2 norms of U and delta U (see "Fast computation of cell similarity",
    // by Paolo Carnevali, dated 5/9/2017).
    const double n = geneCount();
    const double mX = cell.sum1 / n;
    const double sigmaX = sqrt(cell.sum2/n - mX*mX);
    const double zeroX = -mX/sigmaX;    // The X corresponding to a zero x.
    const double mU = cell.sum1LargeExpressionCounts / n;
    const double sigmaU = sqrt(cell.sum2LargeExpressionCounts/n - mU*mU);
    const double zeroU = -mU/sigmaU;    // The U corresponding to a zero u.
    cell.norm2U = (n-storedLargeExpressionCounts.size()) * zeroU*zeroU;    // Contribution of u's which are zero.
    for(const auto& p: storedLargeExpressionCounts) {
        const double& count = p.second;
        const double U = (count-mU) / sigmaU;
        cell.norm2U += U*U;
    }
    cell.norm2U = sqrt(cell.norm2U);
    cell.norm2DeltaU = (n-storedExpressionCounts.size()) * (zeroU-zeroX)*(zeroU-zeroX);    // Contribution of terms in which both x and u are zero.
    // Add the contribution of terms where x and u are not zero, and equal.
    for(const auto& p: storedLargeExpressionCounts) {
        const double& count = p.second;
        const double X = (count-mX) / sigmaX;
        const double U = (count-mU) / sigmaU;
        const double deltaU = X - U;
        cell.norm2DeltaU += deltaU * deltaU;
    }
    // Add the contribution of terms where u is zero but x is not zero.
    for(size_t i=numberToKeep; i<expressionCounts.size(); i++) {
        const auto& p = expressionCounts[i];
        const double& x = p.second;
        const double X = (x-mX) / sigmaX;
        const double U = zeroU;
        const double deltaU = X - U;
        cell.norm2DeltaU += deltaU * deltaU;
    }
    cell.norm2DeltaU = sqrt(cell.norm2DeltaU);


    cell.largestCount = 0.;
    if(!expressionCounts.empty()) {
        cell.largestCount = expressionCounts.front().second;
    }
    cell.largestSmallCount = 0.;
    if(numberToKeep < expressionCounts.size()) {
        cell.largestSmallCount = expressionCounts[numberToKeep].second;
    }
    cell.sum1SmallCounts = 0.;
    for(size_t i=numberToKeep; i<expressionCounts.size(); i++) {
        const auto& p = expressionCounts[i];
        const double& count = p.second;
        cell.sum1SmallCounts += count;
    }

#endif

    // Add this cell to the AllCells set.
    cellSets.cellSets["AllCells"]->push_back(CellId(cells.size()));

    // Store fixed size information for this cell.
    cells.push_back(cell);




    // Sanity checks.
    CZI_ASSERT(cellNames.size() == cells.size());
    CZI_ASSERT(cellMetaData.size() == cells.size());
    CZI_ASSERT(cellExpressionCounts.size() == cells.size());
    CZI_ASSERT(largeCellExpressionCounts.size() == cells.size());
    CZI_ASSERT(cellSets.cellSets["AllCells"]->size() == cells.size());

    // Done.
    return cellNameStringId;
}



// Version of addCell that takes JSON as input.
// The expected JSON can be constructed using Python code modeled from the following:
// import json
// cell = {'metaData': {'CellName': 'abc', 'key1': 'value1'}, 'expressionCounts': {'gene1': 10,'gene2': 20}}
// jSonString = json.dumps(cell)
// expressionMatrix.addCell(json.dumps(jSonString))
// Note the CellName metaData entry is required.
CellId ExpressionMatrix::addCell(
    const string& jsonString,
    size_t maxTermCountForApproximateSimilarityComputation
    )
{

    try {

        // Convert the JSON to a boost::property_tree.
        boost::property_tree::ptree propertyTree;
        std::istringstream jsonStream(jsonString);
        boost::property_tree::read_json(jsonStream, propertyTree);

        // Extract the meta data from the JSON.
        vector< pair<string, string> > metaData;
        const boost::property_tree::ptree& metaDataPropertyTree = propertyTree.get_child("metaData");
        for(const auto& metaDataItem: metaDataPropertyTree) {
            const string& key = metaDataItem.first;
            const boost::property_tree::ptree& valuePropertyTree = metaDataItem.second;
            const string value = valuePropertyTree.get<string>("");
            metaData.push_back(make_pair(key, value));
        }

        // Extract the expression counts from the JSON.
        vector< pair<string, float> > expressionCounts;
        const boost::property_tree::ptree& expressionCountsPropertyTree = propertyTree.get_child("expressionCounts");
        for(const auto& expressionCountsItem: expressionCountsPropertyTree) {
            const string& geneName = expressionCountsItem.first;
            const boost::property_tree::ptree& valuePropertyTree = expressionCountsItem.second;
            const float count = valuePropertyTree.get<float>("");
            expressionCounts.push_back(make_pair(geneName, count));
        }

        // Call the lower level version of addCell.
        return addCell(metaData, expressionCounts, maxTermCountForApproximateSimilarityComputation);

    } catch(...) {

        // If an error occurred, make sure to write out the JSON that caused it.
        cout << "Error processing the following cell JSON:" << endl;
        cout << jsonString << endl;
        throw;
    }

}



// Add cells from data in files with fields separated by commas or by other separators.
// See ExpressionMatrix.hpp for usage information.
// This is new version is more flexible than the old version, ifdef'ed out below.
void ExpressionMatrix::addCells(
    const string& expressionCountsFileName,
    const string& expressionCountsFileSeparators,
    const string& cellMetaDataFileName,
    const string& cellMetaDataFileSeparators,
    size_t maxTermCountForApproximateSimilarityComputation
    )
{
	// Tokenize the cell meta data file and the expression counts file and verify
	// that all lines have the same number of tokens (which cannot be zero), with the possible exception
	// of line one which could have one less token than all other lines.
	vector< vector<string> > cellMetaDataFileLines, expressionCountsFileLines;
	cout << timestamp << "Reading cell meta data file " << cellMetaDataFileName << "." << endl;
	tokenizeFileAndCheck(cellMetaDataFileName, cellMetaDataFileSeparators, cellMetaDataFileLines);
	cout << timestamp << "Reading expression counts file " << expressionCountsFileName << "." << endl;
	tokenizeFileAndCheck(expressionCountsFileName, expressionCountsFileSeparators, expressionCountsFileLines);



	// Get the meta data field names from the header of the cell meta data file.
	// We have to account for the fact that the header line might or might not contain an
	// initial field, which if present is ignored.
	vector<string> metaDataFieldNames;
	{
		size_t skip = 1;
		if(cellMetaDataFileLines[0].size() != cellMetaDataFileLines[1].size()) {
			CZI_ASSERT(cellMetaDataFileLines[0].size() == cellMetaDataFileLines[1].size()-1); // This was checked by tokenizeFileAndCheck.
			skip = 0;
		}
		copy(cellMetaDataFileLines[0].begin()+skip, cellMetaDataFileLines[0].end(), back_inserter(metaDataFieldNames));
	}



	// Check that there are no duplications in the meta data field names.
	{
		set<string> metaDataFieldNamesSet;
		for(const string& metaDataFieldName: metaDataFieldNames) {
			if(metaDataFieldNamesSet.find(metaDataFieldName) != metaDataFieldNamesSet.end()) {
				throw runtime_error("Duplicate meta data field " + metaDataFieldName);
			}
			metaDataFieldNamesSet.insert(metaDataFieldName);
		}
	}



	// Get the cell names from the header of the expression count file, and create a
	// corresponding index map.
	// We have to account for the fact that the header line might or might not contain an
	// initial field, which if present is ignored.
	// Note that not all of these cells will make it into the system:
	// the ones that have no entry in the meta data file will be skipped.
	vector<string> expressionFileCellNames;
	{
		size_t skip = 1;
		if(expressionCountsFileLines[0].size() != expressionCountsFileLines[1].size()) {
			CZI_ASSERT(expressionCountsFileLines[0].size() == expressionCountsFileLines[1].size()-1); // This was checked by tokenizeFileAndCheck.
			skip = 0;
		}
		copy(expressionCountsFileLines[0].begin()+skip, expressionCountsFileLines[0].end(), back_inserter(expressionFileCellNames));
	}
	map<string, CellId> expressionFileCellNamesMap;
	for(size_t i=0; i<expressionFileCellNames.size(); i++) {
		const string& cellName = expressionFileCellNames[i];
		if(expressionFileCellNamesMap.find(cellName) != expressionFileCellNamesMap.end()) {
			throw runtime_error("Cell " + cellName + " has more than one column in the expression counts file.");
		}
		expressionFileCellNamesMap.insert(make_pair(cellName, i));
	}



	// Summarize the number of cells, genes, and meta data names seen in each file.
	cout << "Cell meta data file " << cellMetaDataFileName << " contains data for ";
	cout << cellMetaDataFileLines.size()-1 << " cells and " << metaDataFieldNames.size() << " meta data names." << endl;
	cout << "Expression counts file " << expressionCountsFileName << " contains data for ";
	cout << expressionFileCellNames.size() << " cells and " << expressionCountsFileLines.size()-1 << " genes." << endl;


	// Add the genes.
	// We want to add them independently of the cells, so they all get added, even the ones
	// for which all cells have zero count.
	CZI_ASSERT(expressionCountsFileLines.size() > 0); 	// This was checked by tokenizeFileAndCheck.
	for(size_t i=1; i<expressionCountsFileLines.size(); i++) {
		const vector<string>& line = expressionCountsFileLines[i];
		CZI_ASSERT(line.size() > 0); 	// This was checked by tokenizeFileAndCheck.
		addGene(line.front());
	}



	// Loop over cells in the cell meta data file, but only add the ones that also appear
	// in the expression counts file.
	cout << timestamp << "Storing expression counts and cell meta data." << endl;
	CZI_ASSERT(cellMetaDataFileLines.size() > 1);	// This was checked by tokenizeFileAndCheck.
	CellId addedCellCount = 0;
	for(size_t cellMetaDataFileLine=1; cellMetaDataFileLine<cellMetaDataFileLines.size(); cellMetaDataFileLine++) {
		if((cellMetaDataFileLine%1000) == 0) {
			cout << timestamp << "Working on cell meta data file line " << cellMetaDataFileLine+1 << " of " << cellMetaDataFileLines.size() << endl;
		}
		const vector<string>& metaDataLine = cellMetaDataFileLines[cellMetaDataFileLine];
		CZI_ASSERT(metaDataLine.size() > 1);	// This was checked by tokenizeFileAndCheck.
		const string& cellName = metaDataLine.front();

		// See if this cell appears in the expression counts file.
		// If not, skip this cell.
		const auto it = expressionFileCellNamesMap.find(cellName);
		if(it == expressionFileCellNamesMap.end()) {
			continue;	// It's not in the expression counts file. Skip it.
		}

		// Find the column in the expression counts file that contains data for this cell.
		const size_t expressionFileColumn = it->second + 1;

		// Gather the meta data for this cell.
		CZI_ASSERT(metaDataLine.size() == metaDataFieldNames.size() + 1); // This was checked by tokenizeFileAndCheck.
	    vector< pair<string, string> > thisCellMetaData;
	    thisCellMetaData.push_back(make_pair("CellName", cellName));
	    for(size_t i=0; i<metaDataFieldNames.size(); i++) {
	    	thisCellMetaData.push_back(make_pair(metaDataFieldNames[i], metaDataLine[i+1]));
	    }

	    // Gather the expression counts for this cell.
	    vector< pair<string, float> > thisCellExpressionCounts;
		for(size_t i=1; i<expressionCountsFileLines.size(); i++) {
			const vector<string>& line = expressionCountsFileLines[i];
			CZI_ASSERT(line.size() > 0); 	// This was checked by tokenizeFileAndCheck.
			const string& geneName = line.front();
			const string& expressionCountString = line[expressionFileColumn];
			float expressionCount;
			try {
				expressionCount = lexical_cast<float>(expressionCountString);
			} catch(boost::bad_lexical_cast) {
				throw runtime_error("Invalid expression count " + expressionCountString +
					" for cell " + cellName + " gene " + geneName);
			}
			if(expressionCount != 0.) {
				thisCellExpressionCounts.push_back(make_pair(geneName, expressionCount));
			}
		}

		// Now we can add this cell.
		++addedCellCount;
		addCell(thisCellMetaData, thisCellExpressionCounts, maxTermCountForApproximateSimilarityComputation);
	}

	cout << timestamp << "Added " << addedCellCount << " cells that appear in both the cell meta data file and the expression counts file." << endl;
	cout << "There are " << cellCount() << " cells and " << geneCount() << " genes." << endl;

}



#if 0
// Add cells from data in files with fields separated by commas or by other separators.
// See ExpressionMatrix.hpp for usage information.
// This is the old version that requires the expression count file and the meta data file
// to have exactly the same cells, in the same order.
void ExpressionMatrix::addCells(
    const string& expressionCountsFileName,
    const string& expressionCountsFileSeparators,
    const string& metaDataFileName,
    const string& metaDataFileSeparators,
    size_t maxTermCountForApproximateSimilarityComputation
    )
{

    // Variables used to hold a line of an input file, and its tokenized version
    string line;
    vector<string> tokens;

    // Open the expression count file.
    ifstream expressionCountsFile(expressionCountsFileName);
    if(!expressionCountsFile) {
        throw runtime_error("Error opening the expression count file " + expressionCountsFileName);
    }

    // Get the cell names from the first row of the expression count file.
    getline(expressionCountsFile, line);
    if(!expressionCountsFile) {
        throw runtime_error("Error reading the first line of the expression count file.");
    }
    if(line.empty()) {
        throw runtime_error("The first line of the expression count file is empty.");
    }
    tokenize(expressionCountsFileSeparators, line, tokens);
    if(tokens.size() < 2) {
        throw runtime_error("The first line of the expression count file does not contain the specified separator.");
    }
    vector<string> cellNames(tokens.begin()+1, tokens.end());



    // Get the gene names and the expression counts from the rest of the
    // expression count file.
    vector<string> geneNames;
    vector< vector<float> > counts; // Here, one vector for each gene.
    while(true) {
        if(!counts.empty() && ((counts.size()%1000)==0)) {
            cout << timestamp << "Read expression counts for " << counts.size() << " genes." << endl;
        }
        getline(expressionCountsFile, line);
        if(!expressionCountsFile) {
            break;
        }
        tokenize(expressionCountsFileSeparators, line, tokens);
        if(tokens.size() != cellNames.size()+1) {
            cout << "Unexpected number of tokens in expression counts line." << endl;
            cout << "Expected " << cellNames.size()+1 << " tokens." << endl;
            cout << "Found " << tokens.size() << " tokens." << endl;
            cout << "Offending line:" << endl;
            cout << line << endl;
            cout << "Unexpected number of tokens in expression counts line." << endl;
            cout << "Expected " << cellNames.size()+1 << " tokens." << endl;
            cout << "Found " << tokens.size() << " tokens." << endl;
            if(tokens.size() == cellNames.size()+2) {
            	cout << "It is possible that the header line of the expression file "
            		"is missing the first (ignored) token above the column containing gene names." << endl;
            }
            throw runtime_error("Unexpected number of tokens in expression counts line.");
        }
        const string& geneName = tokens.front();
        geneNames.push_back(geneName);
        addGene(geneName);
        counts.resize(counts.size()+1, vector<float>(cellNames.size()));
        try {
            for(size_t i=0; i<cellNames.size(); i++) {
                counts.back()[i] = lexical_cast<float>(tokens[i+1]);
            }
        } catch(bad_lexical_cast) {
            cout << "Error extracting expression counts from expression count line:" << endl;
            cout << line << endl;
            throw runtime_error("Error extracting expression counts from expression count line.");
        }
    }



    // Vectors to contain meta data names (the same for all cells)
    // and values (different from each cell).
    // Initialize them with just the cell names.
    vector<string> metaDataNames;
    vector< vector<string> > metaDataValues(cellNames.size());
    metaDataNames.push_back("CellName");
    for(CellId cellId=0; cellId<cellNames.size(); cellId++) {
        metaDataValues[cellId].push_back(cellNames[cellId]);
    }



    // If an input meta data file was specified, read cell meta data from it.
    if(!metaDataFileName.empty()) {

        // Open the meta data file.
        ifstream metaDataFile(metaDataFileName);
        if(!metaDataFile) {
            throw runtime_error("Error opening cell meta data file " + metaDataFileName);
        }

        // Read the meta data names, which are the same for all cells,
        // from the first line of the meta data file.
        getline(metaDataFile, line);
        if(!metaDataFile) {
            throw runtime_error("Error reading first row of meta data file.");
        }
        tokenize(metaDataFileSeparators, line, tokens);
        if(tokens.size() < 2) {
            throw runtime_error("Unexpected format of first line of meta data file. It is possible that the incorrect separator was specified.");
        }
        metaDataNames.insert(metaDataNames.end(), tokens.begin()+1, tokens.end());

        // Verify that there are no duplications in the meta data names.
        for(size_t i=1; i<metaDataNames.size(); i++) {
            const string& iName = metaDataNames[i];
            for(size_t j=0; j<i; j++) {
                if(metaDataNames[j] == iName) {
                    throw runtime_error("Duplicate meta data name " + iName);
                }
            }
        }

        // Read cell meta data from the rest of the cell meta data file.
        for(CellId cellId=0; cellId<cellNames.size(); cellId++) {
            getline(metaDataFile, line);
            if(!metaDataFile) {
                throw runtime_error("Error reading meta data file line for cell " + cellNames[cellId]);
            }
            tokenize(metaDataFileSeparators, line, tokens);
            if(tokens.size() != metaDataNames.size()) { // metaDataNames also contains CellName.
                cout << "Unexpected number of tokens in meta data file line:" << endl;
                cout << line << endl;
                cout << "Expected " << metaDataNames.size() << " tokens ";
                cout << " for " << metaDataNames.size()-1 << " meta data items, but got ";
                cout << tokens.size()  << " tokens." << endl;
                cout << "Meta data names:" << endl;
                for(size_t i=1; i<metaDataNames.size(); i++) {
                    cout << i << " " << metaDataNames[i] << endl;
                }
                cout << "Unexpected number of tokens in meta data file line:" << endl;
                cout << line << endl;
                // Write this portion of the message again, for better visibility.
                cout << "Expected " << metaDataNames.size() << " tokens ";
                cout << " for " << metaDataNames.size()-1 << " meta data items, but got ";
                cout << tokens.size()  << " tokens." << endl;
                if(tokens.size() == metaDataNames.size()+1) {
                	cout << "It is possible that the meta data file is missing "
                		"the first (ignored) field of the first line (the field above the cell name)." << endl;
                }
                throw runtime_error("Unexpected number of tokens in meta data file line.");
            }
            if(tokens.front() != cellNames[cellId]) {
                cout << "Expected the following cell name in line of cell data file: " << cellNames[cellId];
                cout << " but found " << tokens.front() << "." << endl;
                cout << " Offending line is:" << endl;
                cout << line << endl;
                throw runtime_error("Unexpected cell name in meta data file.");
            }
            metaDataValues[cellId].insert(metaDataValues[cellId].end(), tokens.begin()+1, tokens.end());
        }
    }


    // Count the number of genes that have zero counts for all cells.
    GeneId zeroCount = 0;
    for(GeneId geneId=0; geneId<geneNames.size(); geneId++) {
        bool foundNonZero = false;
        for(CellId cellId=0; cellId!=cellNames.size(); cellId++) {
            if(counts[geneId][cellId] != 0) {
                foundNonZero = true;
                break;
            }
        }
        if(!foundNonZero) {
            ++zeroCount;
        }
    }
    cout << "Found " << zeroCount << " genes with zero counts for all cells." << endl;



    // Now we have all the information we need to add the cells one by one.
    vector< pair<string, string> > thisCellMetaData;
    vector< pair<string, float> > thisCellExpressionCounts;
    for(CellId cellId=0; cellId!=cellNames.size(); cellId++) {

        // Fill in the meta data.
        CZI_ASSERT(metaDataValues[cellId].size() == metaDataNames.size());
        for(size_t i=0; i<metaDataNames.size(); i++) {
            thisCellMetaData.push_back(make_pair(metaDataNames[i], metaDataValues[cellId][i]));
        }

        // Fill in the expression counts.
        for(GeneId geneId=0; geneId<geneNames.size(); geneId++) {
            const float thisCount = counts[geneId][cellId];
            if(thisCount != 0.) {
                thisCellExpressionCounts.push_back(make_pair(geneNames[geneId], thisCount));
            }
        }

        addCell(thisCellMetaData, thisCellExpressionCounts, maxTermCountForApproximateSimilarityComputation);
        thisCellMetaData.clear();
        thisCellExpressionCounts.clear();
    }


    cout << "The expression matrix has " << geneCount();
    cout << " genes and " << cellCount() << " cells." << endl;
    cout << "The total number of expression counts is " << cellExpressionCounts.totalSize() << endl;
    cout << "The total number of large expression counts is " << largeCellExpressionCounts.totalSize() << endl;
}
#endif



// Return a cell id given a string.
// The string can be a cell name or a CellId (an integer).
// Returns invalidCellId if the cell was not found.
CellId ExpressionMatrix::cellIdFromString(const string& s)
{
    // If the strings represent a CellId in the expected range, treat it as a cell id.
    try {
        const CellId cellId = lexical_cast<CellId>(s);
        if(cellId < cellCount()) {
            return cellId;
        }
    } catch(bad_lexical_cast) {
        // Nothing to worry about. The string was not a number.
    }

    // Not an integer. Treat it as a cell name.
    return cellNames(s);
}



// Return a gene id given a string.
// The string can be a gene name or GeneId (a string).
// Returns imnvalidGeneId if the gene was not found.
GeneId ExpressionMatrix::geneIdFromString(const string& s)
{
    // If the strings represent a GeneId in the expected range, treat it as a gene id.
    try {
        const GeneId geneId = lexical_cast<GeneId>(s);
        if(geneId < geneCount()) {
            return geneId;
        }
    } catch(bad_lexical_cast) {
        // Nothing to worry about. The string was not a number.
    }

    // Not an integer. Treat it as a gene name.
    return geneNames(s);
}



// Return the value of a specified meta data field for a given cell.
// Returns an empty string if the cell does not have the specified meta data field.
string ExpressionMatrix::getMetaData(CellId cellId, const string& name) const
{
    // Find the string id of the name.
    // If it does not exist, return an empty string.
    const StringId nameId = cellMetaDataNames(name);
    if(nameId == cellMetaDataNames.invalidStringId) {
        return "";
    }
    return getMetaData(cellId, nameId);
}
string ExpressionMatrix::getMetaData(CellId cellId, StringId nameId) const
{

    // Scan the name/value pairs for this cell, looking for nameId.
    for(const auto& metaDataPair: cellMetaData[cellId]) {
        if(metaDataPair.first == nameId) {
            const StringId valueId = metaDataPair.second;
            if(valueId == cellMetaDataValues.invalidStringId) {
                return "";  // Should never happen, but just in case.
            } else {
                return cellMetaDataValues[valueId];
            }
        }
    }

    // We did not find it. Return an empty string.
    return "";
}



// Set a meta data (name, value) pair for a given cell.
// If the name already exists for that cell, the value is replaced.
void ExpressionMatrix::setMetaData(CellId cellId, const string& name, const string& value)
{
    const StringId nameId = cellMetaDataNames[name];
    const StringId valueId = cellMetaDataValues[value];
    setMetaData(cellId, nameId, valueId);
}
void ExpressionMatrix::setMetaData(CellId cellId, StringId nameId, const string& value)
{
    const StringId valueId = cellMetaDataValues[value];
    setMetaData(cellId, nameId, valueId);
}
void ExpressionMatrix::setMetaData(CellId cellId, StringId nameId, StringId valueId)
{

    // Scan the existing meta data for this cell, looking for this name.
    for(auto& p: cellMetaData[cellId]) {
        if(p.first == nameId) {
            p.second = valueId; // The name already exists. replace the value.
            return;
        }
    }

    // The name did not exist for this cell. Add this (name, value) pair.
    cellMetaData.push_back(cellId, make_pair(nameId, valueId));
    incrementCellMetaDataNameUsageCount(nameId);
}



void ExpressionMatrix::incrementCellMetaDataNameUsageCount(StringId nameId)
{
    if(cellMetaDataNamesUsageCount.size() <= nameId) {
        // This is a new name.
        // cout << "***A " << cellMetaDataNamesUsageCount.size() << " " << nameId+1);
        CZI_ASSERT(cellMetaDataNamesUsageCount.size() == nameId);
        cellMetaDataNamesUsageCount.push_back(1);
    } else {

        // This is an existing name.
        ++(cellMetaDataNamesUsageCount[nameId]);
    }

}



void ExpressionMatrix::decrementCellMetaDataNameUsageCount(StringId nameId)
{
    CZI_ASSERT(nameId < cellMetaDataNamesUsageCount.size());
    CZI_ASSERT(cellMetaDataNamesUsageCount[nameId] > 0);
    --(cellMetaDataNamesUsageCount[nameId]);
}



// Return the raw expression count for a given CellId and GeneId.
float ExpressionMatrix::getExpressionCount(CellId cellId, GeneId geneId) const
{
    CZI_ASSERT(cellId < cellCount());
    CZI_ASSERT(geneId < geneCount());
    const auto& counts = cellExpressionCounts[cellId];
    auto it = lower_bound(counts.begin(), counts.end(), pair<int, float>(geneId, 0.), OrderPairsByFirstOnly< pair<int, float> >());
    if(it == counts.end() || it->first != geneId) {
        return 0.;
    } else {
        return it->second;
    }
}



// Compute the average expression vector for a given set of cells.
// The last parameter controls the normalization used for the expression counts
// for averaging:
// 0: no normalization (raw read counts).
// 1: L1 normalization (fractional read counts).
// 2: L2 normalization.
void ExpressionMatrix::computeAverageExpression(
	const vector<CellId> cellIds,
	vector<double>& averageExpression,
	NormalizationMethod normalizationMethod) const
{
	// Initialize the average expression to zero.
	averageExpression.resize(geneCount());
	fill(averageExpression.begin(), averageExpression.end(), 0.);



	// Accumulate the contribution of all the cells.
	for(const CellId cellId: cellIds) {

		// Compute the normalization factor for this cell.
		const Cell& cell = cells[cellId];
		double factor;
		switch(normalizationMethod) {
		case NormalizationMethod::None:
			factor = 1.;
			break;
		case NormalizationMethod::L1:
			factor = cell.norm1Inverse;
			break;
		case NormalizationMethod::L2:
			factor = cell.norm2Inverse;
			break;
		default:
			CZI_ASSERT(0);
		}


		// Add all of the expression counts for this cell.
		for(const auto& p: cellExpressionCounts[cellId]) {
			const GeneId geneId = p.first;
			float count = p.second;
			const double normalizedCount = factor * double(count);
			averageExpression[geneId] += normalizedCount;
		}
	}



	// Divide by the number of cells.
	const double factor = 1./double(cellIds.size());
	for(double& a: averageExpression) {
		a *= factor;
	}



	// Normalize as requested.
	switch(normalizationMethod) {
	case NormalizationMethod::None:
		break;
	case NormalizationMethod::L1:
	{
		const double factor = 1. / std::accumulate(averageExpression.begin(), averageExpression.end(), 0.);
		for(double& a: averageExpression) {
			a *= factor;
		}
		break;
	}
	case NormalizationMethod::L2:
	{
		double sum = 0.;
		for(const double& a: averageExpression) {
			sum += a*a;;
		}
		const double factor = 1./sqrt(sum);
		for(double& a: averageExpression) {
			a *= factor;
		}
		break;
	}
	default:
		CZI_ASSERT(0);
	}

}



// Compute a sorted histogram of a given meta data field.
void ExpressionMatrix::histogramMetaData(
    const CellSets::CellSet& cellSet,
    StringId metaDataNameId,
    vector< pair<string, size_t> >& sortedHistogram) const
{
    // Create the histogram.
    map<string, size_t> histogram;
    for(const CellId cellId: cellSet) {
        const string metaDataValue = getMetaData(cellId, metaDataNameId);
        const auto it = histogram.find(metaDataValue);
        if(it == histogram.end()) {
            histogram.insert(make_pair(metaDataValue, 1));
        } else {
            ++(it->second);
        }
    }


    // Sort the histogram by decreasing frequency.
    sortedHistogram.clear();
    copy(histogram.begin(), histogram.end(), back_inserter(sortedHistogram));
    sort(sortedHistogram.begin(), sortedHistogram.end(), OrderPairsBySecondGreaterThenByFirstLess< pair<string, size_t> >());
}



// Compute the similarity between two cells given their CellId.
// The similarity is the correlation coefficient of their
// expression counts.
double ExpressionMatrix::computeCellSimilarity(CellId cellId0, CellId cellId1) const
{
    // Compute the scalar product of the expression counts for the two cells.
    typedef pair<GeneId, float>const* Iterator;
    const Iterator begin0 = cellExpressionCounts.begin(cellId0);
    const Iterator end0 = cellExpressionCounts.end(cellId0);
    const Iterator begin1 = cellExpressionCounts.begin(cellId1);
    const Iterator end1 = cellExpressionCounts.end(cellId1);
    Iterator it0 = begin0;
    Iterator it1 = begin1;
    double scalarProduct = 0.;
    while((it0 != end0) && (it1 != end1)) {
        const GeneId geneId0 = it0->first;
        const GeneId geneId1 = it1->first;

        if(geneId0 < geneId1) {
            ++it0;
        } else if(geneId1 < geneId0) {
            ++it1;
        } else {
            scalarProduct += it0->second * it1->second;
            ++it0;
            ++it1;
        }
    }

    // Compute the correlation coefficient.
    // See, for example, https://en.wikipedia.org/wiki/Correlation_and_dependence
    const double n = geneCount();
    const Cell& cell0 = cells[cellId0];
    const Cell& cell1 = cells[cellId1];
    const double numerator = n*scalarProduct - cell0.sum1*cell1.sum1;
    const double denominator = sqrt(
            (n*cell0.sum2 - cell0.sum1*cell0.sum1) *
            (n*cell1.sum2 - cell1.sum1*cell1.sum1)
            );
    return numerator / denominator;
}



// Approximate but fast computation of the similarity between two cells.
double
    ExpressionMatrix::computeApproximateCellSimilarity(CellId cellId0, CellId cellId1) const
{
    // The lower bound is computed just like the exact similarity, but we only use the
    // largeCellExpressionCounts instead of the cellExpressionCounts.

    // Compute the scalar product of the large expression counts for the two cells.
    typedef pair<GeneId, float>const* Iterator;
    const Iterator begin0 = largeCellExpressionCounts.begin(cellId0);
    const Iterator end0 = largeCellExpressionCounts.end(cellId0);
    const Iterator begin1 = largeCellExpressionCounts.begin(cellId1);
    const Iterator end1 = largeCellExpressionCounts.end(cellId1);
    Iterator it0 = begin0;
    Iterator it1 = begin1;
    double scalarProduct = 0.;
    while((it0 != end0) && (it1 != end1)) {
        const GeneId geneId0 = it0->first;
        const GeneId geneId1 = it1->first;

        if(geneId0 < geneId1) {
            ++it0;
        } else if(geneId1 < geneId0) {
            ++it1;
        } else {
            scalarProduct += it0->second * it1->second;
            ++it0;
            ++it1;
        }
    }


    // Compute the correlation coefficient.
    // See, for example, https://en.wikipedia.org/wiki/Correlation_and_dependence
    const double n = geneCount();
    const Cell& cell0 = cells[cellId0];
    const Cell& cell1 = cells[cellId1];
    const double numerator = n*scalarProduct - cell0.sum1LargeExpressionCounts*cell1.sum1LargeExpressionCounts;
    const double denominator = sqrt(
            (n*cell0.sum2LargeExpressionCounts - cell0.sum1LargeExpressionCounts*cell0.sum1LargeExpressionCounts) *
            (n*cell1.sum2LargeExpressionCounts - cell1.sum1LargeExpressionCounts*cell1.sum1LargeExpressionCounts)
            );
    return numerator / denominator;


#if 0
    // Additional code used to compute bounds on approximate cell similarity.
    // None of the bounds turned out to be sufficiently tight, so we are not using these.

    // The lower bound is computed using the exact expression, replacing the scalar product
    // of all terms with the scalar product of the large terms.
    // The upper bound is computed in the same way, using the upper bound for the scalar product.
    {
        const double numerator = n*scalarProduct - cell0.sum1*cell1.sum1;
        const double denominator = sqrt(
                (n*cell0.sum2 - cell0.sum1*cell0.sum1) *
                (n*cell1.sum2 - cell1.sum1*cell1.sum1)
                );
        approximateCellSimilarity.lowerBound1 = numerator / denominator;
        const double scalarProductErrorBound =
            cell0.largeExpressionCountsNorm2 * cell1.smallExpressionCountsNorm2 +
            cell1.largeExpressionCountsNorm2 * cell0.smallExpressionCountsNorm2 +
            cell0.smallExpressionCountsNorm2 * cell0.smallExpressionCountsNorm2;
        approximateCellSimilarity.upperBound1 = (numerator + n*scalarProductErrorBound) / denominator;
    }


    const double errorEstimate = (cell0.norm2U*cell1.norm2DeltaU + cell1.norm2U*cell0.norm2DeltaU + cell0.norm2DeltaU*cell1.norm2DeltaU) / n;
    approximateCellSimilarity.lowerBound2 = approximateCellSimilarity.estimate - errorEstimate;
    approximateCellSimilarity.upperBound2 = approximateCellSimilarity.estimate + errorEstimate;

    cout << "Cell0: " << cell0.norm2U << " " << cell0.norm2DeltaU << endl;
    cout << "Cell1: " << cell1.norm2U << " " << cell1.norm2DeltaU << endl;


    // Upper bound on page 3 of
    // "Fast computation of cell similarity", by Paolo Carnevali, dated 5/8/2017.
    {
        const double numerator = n*scalarProduct - cell0.sum1*cell1.sum1;
        const double denominator = sqrt(
            (n*cell0.sum2 - cell0.sum1*cell0.sum1) *
            (n*cell1.sum2 - cell1.sum1*cell1.sum1)
            );
        const double scalarProductErrorBound =
            cell0.largestCount * cell1.sum1SmallCounts +
            cell1.largestCount * cell0.sum1SmallCounts +
            min(cell0.largestSmallCount * cell1.sum1SmallCounts, cell1.largestSmallCount * cell0.sum1SmallCounts);
        approximateCellSimilarity.upperBound3 = (numerator + n*scalarProductErrorBound) / denominator;
    }


    return approximateCellSimilarity;
#endif
}



// Create a new gene set consisting of genes whose name matches a given reguklar expression.
bool ExpressionMatrix::createGeneSetUsingGeneNames(const string& geneSetName, const string& regexString)
{
	// Check if a gene set with this name already exists.
	if(geneSets.find(geneSetName) != geneSets.end()) {
		return false;
	}

	// Create the regular expression we are going to match.
    const boost::regex regex(regexString);

	// Create the new gene set.
	GeneSet& geneSet = geneSets[geneSetName];
	geneSet.createNew(directoryName + "/GeneSet-" + geneSetName);
	for(GeneId geneId=0; geneId!=geneCount(); geneId++) {
		const string geneName = geneNames[geneId];
		if(boost::regex_match(geneName, regex)) {
			geneSet.addGene(geneId);
		}
	}


	return true;
}



bool ExpressionMatrix::createGeneSetIntersection(const string& inputSetsNames, const string& outputSetName)
{
    return createGeneSetIntersectionOrUnion(inputSetsNames, outputSetName, false);
}
bool ExpressionMatrix::createGeneSetUnion(const string& inputSetsNames, const string& outputSetName)
{
    return createGeneSetIntersectionOrUnion(inputSetsNames, outputSetName, true);
}
bool ExpressionMatrix::createGeneSetIntersectionOrUnion(
	const string& commaSeparatedInputSetsNames,
	const string& outputSetName,
	bool doUnion)
{
    // See if a gene set with the name of the output gene set already exists.
    if(geneSets.find(outputSetName) != geneSets.end()) {
        cout << "Gene set " << outputSetName << " already exists." << endl;
        return false;
    }

    // Parse the input gene sets.
    vector<string> inputSetsNames;
    boost::algorithm::split(inputSetsNames, commaSeparatedInputSetsNames, boost::is_any_of(","));

    // Check that all input gene sets exist.
    for(const string& inputSetName: inputSetsNames) {
        if(geneSets.find(inputSetName) == geneSets.end()) {
            cout << "gene set " << inputSetName << " does not exists." << endl;
            return false;
        }
    }

    // Compute the intersection or union.
    vector<GeneId> outputSetGenes;
    for(size_t i=0; i<inputSetsNames.size(); i++) {
        const string& inputSetName = inputSetsNames[i];
        vector<GeneId> inputSetGenes;
        geneSets[inputSetName].getSortedGenes(inputSetGenes);
        if(i == 0) {
            outputSetGenes = inputSetGenes;
        } else {
            vector<GeneId> newOutputSetGenes;
            if(doUnion) {
                std::set_union(
                    outputSetGenes.begin(), outputSetGenes.end(),
                    inputSetGenes.begin(), inputSetGenes.end(),
                    back_inserter(newOutputSetGenes));
            } else {
                std::set_intersection(
                    outputSetGenes.begin(), outputSetGenes.end(),
                    inputSetGenes.begin(), inputSetGenes.end(),
                    back_inserter(newOutputSetGenes));
            }
            outputSetGenes.swap(newOutputSetGenes);
        }
    }



    // Store this gene set.
    GeneSet& outputGeneSet = geneSets[outputSetName];
    outputGeneSet.createNew(directoryName + "/GeneSet-" + outputSetName);
    for(const GeneId geneId: outputSetGenes) {
        outputGeneSet.addGene(geneId);
    }

    return true;
}



// Create a new cell set that contains cells for which
// the value of a specified meta data field matches
// a given regular expression.
// Return true if successful.
bool ExpressionMatrix::createCellSetUsingMetaData(
    const string& cellSetName,          // The name of the cell set to be created.
    const string& metaDataFieldName,    // The name of the meta data field to be used.
    const string& regexString           // The regular expression that must be matched for a cell to be added to the set.
    )
{
    // See if a cell set with this name already exists.
    if(cellSets.exists(cellSetName)) {
        cout << "Cell set " << cellSetName << " already exists." << endl;
        return false;
    }

    // Create the regular expression we are going to match.
    const boost::regex regex(regexString);



    // Find the cells that belong to the new cell set.
    vector<CellId> cellSet;

    // Loop over all cells.
    for(CellId cellId=0; cellId<cells.size(); cellId++) {

        // Loop over the meta data fields for this cell.
        const auto metaData = cellMetaData[cellId];
        for(const pair<StringId, StringId>& p: metaData) {

            // If the meta data name is not the specified one, skip.
            if(!cellMetaDataNames.equal(p.first, metaDataFieldName)) {
                continue;
            }

            // If the meta data value matches the regular expression,
            // add the cell to the set and stop looping over the remaining
            // meta data fields.
            const auto metaDataValue = cellMetaDataValues(p.second);
            if(boost::regex_match(metaDataValue.begin(), metaDataValue.end(), regex)) {
                cellSet.push_back(cellId);
                break;
            }
        }
    }



    // Store this cell set.
    cellSets.addCellSet(cellSetName, cellSet);
    // cout << "New cell set " << cellSetName << " contains " << cellSet.size() << " cells." << endl;

    return true;
}




// Create a new cell set as the intersection or union of two or more existing cell sets.
// The input cell sets are specified comma separated in the first argument.
// Return true if successful, false if one of the input cell sets does not exist
// or the output cell set already exists.
// All sets are stored sorted.
bool ExpressionMatrix::createCellSetIntersection(const string& inputSetsNames, const string& outputSetName)
{
    return createCellSetIntersectionOrUnion(inputSetsNames, outputSetName, false);
}
bool ExpressionMatrix::createCellSetUnion(const string& inputSetsNames, const string& outputSetName)
{
    return createCellSetIntersectionOrUnion(inputSetsNames, outputSetName, true);
}
bool ExpressionMatrix::createCellSetIntersectionOrUnion(const string& commaSeparatedInputSetsNames, const string& outputSetName, bool doUnion)
{
    // See if a cell set with the name of the output cell set already exists.
    if(cellSets.exists(outputSetName)) {
        cout << "Cell set " << outputSetName << " already exists." << endl;
        return false;
    }

    // Parse the input cell sets.
    vector<string> inputSetsNames;
    boost::algorithm::split(inputSetsNames, commaSeparatedInputSetsNames, boost::is_any_of(","));

    // Check that all input cell sets exist.
    for(const string& inputSetName: inputSetsNames) {
        if(!cellSets.exists(inputSetName)) {
            cout << "Cell set " << inputSetName << " does not exists." << endl;
            return false;
        }
    }

    // Compute the intersection or union.
    vector<CellId> outputSet;
    for(size_t i=0; i<inputSetsNames.size(); i++) {
        const string& inputSetName = inputSetsNames[i];
        const auto& inputSet = *cellSets.cellSets[inputSetName];
        if(i == 0) {
            outputSet.insert(outputSet.end(), inputSet.begin(), inputSet.end());
        } else {
            vector<CellId> newOutputSet;
            if(doUnion) {
                std::set_union(
                    outputSet.begin(), outputSet.end(),
                    inputSet.begin(), inputSet.end(),
                    back_inserter(newOutputSet));
            } else {
                std::set_intersection(
                    outputSet.begin(), outputSet.end(),
                    inputSet.begin(), inputSet.end(),
                    back_inserter(newOutputSet));
            }
            outputSet.swap(newOutputSet);
        }
    }



    // Store this cell set.
    cellSets.addCellSet(outputSetName, outputSet);
    // cout << "New cell set " << outputSetName << " contains " << outputSet.size() << " cells." << endl;

    return true;
}



bool ExpressionMatrix::createCellSetDifference(
	const string& inputSetName0,
	const string& inputSetName1,
	const string& outputSetName)
{
    // See if a cell set with the name of the output cell set already exists.
    if(cellSets.exists(outputSetName)) {
        cout << "Cell set " << outputSetName << " already exists." << endl;
        return false;
    }



    // Locate the input cell sets.
    const auto it0 = cellSets.cellSets.find(inputSetName0);
    if(it0 == cellSets.cellSets.end()) {
        cout << "Cell set " << inputSetName0 << " does not exists." << endl;
    	return false;
    }
    const CellSets::CellSet& inputSet0 = *(it0->second);
    const auto it1 = cellSets.cellSets.find(inputSetName1);
    if(it1 == cellSets.cellSets.end()) {
        cout << "Cell set " << inputSetName1 << " does not exists." << endl;
    	return false;
    }
    const CellSets::CellSet& inputSet1 = *(it1->second);



    // Compute the difference.
    vector<CellId> outputSet;
	std::set_difference(
		inputSet0.begin(), inputSet0.end(),
		inputSet1.begin(), inputSet1.end(),
		back_inserter(outputSet));



    // Store this cell set.
    cellSets.addCellSet(outputSetName, outputSet);
    return true;
}



// Create a new cell set by downsampling an existing cell set.
bool ExpressionMatrix::downsampleCellSet(
	const string& inputCellSetName,
	const string& outputCellSetName,
	double probability,
	int seed)
{

	// Locate the input cell set.
	const auto it = cellSets.cellSets.find(inputCellSetName);
	if(it == cellSets.cellSets.end()) {
		return false;
	}
    const CellSets::CellSet& inputCellSet = *(it->second);

    // Create the new cell set.
    vector<CellId> outputCellSet;

	// Prepare to generate uniformly distributed numbers between 0 and 1.
	using RandomSource = boost::mt19937;
	using UniformDistribution = boost::uniform_01<>;
	RandomSource randomSource(seed);
	UniformDistribution uniformDistribution;
	boost::variate_generator<RandomSource, UniformDistribution> uniformGenerator(randomSource, uniformDistribution);


    // Loop over all cells in the input cell set.
    // Add each one of them to the output cell set with the specified probability.
    for(const CellId cellId: inputCellSet) {
    	if(uniformGenerator() < probability) {
    		outputCellSet.push_back(cellId);
    	}
    }



    // Store the new cell set.
    cellSets.addCellSet(outputCellSetName, outputCellSet);

    return true;
}



// Create a new graph.
// Graphs are not persistent (they are stored in memory only).
void ExpressionMatrix::createCellSimilarityGraph(
    const string& graphName,            // The name of the graph to be created. This is used as a key in the graph map.
    const string& cellSetName,          // The cell set to be used.
    const string& similarPairsName,     // The name of the SimilarPairs object to be used to create the graph.
    double similarityThreshold,         // The minimum similarity to create an edge.
    size_t maxConnectivity              // The maximum number of neighbors (k of the k-NN graph).
 )
{
    // A graph with this name should not already exist.
    if(graphs.find(graphName) != graphs.end()) {
        throw runtime_error("Graph " + graphName + " already exists.");
    }

    // Locate the cell set.
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exists.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);

    // Create the graph.
    typedef boost::shared_ptr<CellSimilarityGraph> GraphSharedPointer;
    const GraphSharedPointer graph = GraphSharedPointer(new CellSimilarityGraph(
        cellSet,
        directoryName + "/SimilarPairs-" + similarPairsName,
        similarityThreshold,
        maxConnectivity
        ));

    // Create the GraphInformation object that will be stored with the graph.
    GraphInformation graphInformation;
    graphInformation.cellSetName = cellSetName;
	graphInformation.similarPairsName = similarPairsName;
	graphInformation.similarityThreshold = similarityThreshold;
	graphInformation.maxConnectivity = maxConnectivity;

	// Remove isolated vertices.
	graphInformation.isolatedVertexCount = graph->removeIsolatedVertices();
	graphInformation.vertexCount = num_vertices(*graph);
	graphInformation.edgeCount = num_edges(*graph);

    // Store it.
    graphs.insert(make_pair(graphName, make_pair(graphInformation, graph)));

}



// Store the cluster ids in a graph in a meta data field.
void ExpressionMatrix::storeClusterId(
    const string& metaDataName,
    const CellSimilarityGraph& graph)
{
    // Find the string id corresponding to the specified meta data name.
    // This adds it to the table if not already present.
    const StringId metaDataNameStringId = cellMetaDataNames[metaDataName];

    // Loop over all vertices in the graph.
    BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
        const CellSimilarityGraphVertex& vertex = graph[v];

        // Extract the cell id and the cluster id.
        const CellId cellId = vertex.cellId;
        const uint32_t clusterId = vertex.clusterId;

        // Store the cluster id as cell meta data.
        // If the name already exists for this cell, the value is replaced.
        setMetaData(cellId, metaDataNameStringId, lexical_cast<string>(clusterId));
    }
}



// Compute gene information content in bits for a given gene set and cell set,
// using the specified normalization method.
// We do it one gene at a time to avoid the need for an amount of
// memory proportional to the product of the number of cells
// times the number of genes.
void ExpressionMatrix::computeGeneInformationContent(
	const GeneSet& geneSet,
	const CellSets::CellSet& cellSet,
	NormalizationMethod normalizationMethod,
	vector<float>& geneInformationContent) const
{
	geneInformationContent.reserve(geneSet.size());
	geneInformationContent.clear();
	for(const GeneId geneId: geneSet) {
		if(geneId>0 && (geneId%1000)==0) {
			cout << geneId << endl;
		}
		geneInformationContent.push_back(computeGeneInformationContent(geneId, cellSet, normalizationMethod));
	}

}


float ExpressionMatrix::computeGeneInformationContent(
	GeneId geneId,
	const CellSets::CellSet& cellSet,
	NormalizationMethod normalizationMethod) const
{

	// Create a vector of expression counts for this gene and for all cells in the cell set,
	// using the requested normalization.
	// Note that we use the normalization defined using all genes.
	vector<float> count;
	count.reserve(cellSet.size());
	for(const CellId cellId: cellSet) {
		const Cell& cell = cells[cellId];
		float c = getExpressionCount(cellId, geneId);
		switch(normalizationMethod) {
		case NormalizationMethod::L1:
			c *= float(cell.norm1Inverse);
			break;
		case NormalizationMethod::L2:
			c *= float(cell.norm2Inverse);
			break;
		default:
			break;
		}
		count.push_back(c);
	}



	// Compute the sum, using double precision.
	double sum = 0.;
	for(const float c: count) {
		sum += double(c);
	}

	// Compute the information content.
	double informationContent = log(double(cellSet.size())); // Equally distributed.
	const double inverseSum = 1./sum;	// No problem with division by zero - never used if sum is zero
	for(const float c: count) {
		if(c > 0.) {
			const double p = c * inverseSum;
			informationContent += p*log(p);
		}
	}


	// Convert to bits.
	informationContent /= log(2.);

	return float(informationContent);
}
