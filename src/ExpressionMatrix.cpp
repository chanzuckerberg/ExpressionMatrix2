#include "ExpressionMatrix.hpp"
#include "CellGraph.hpp"
#include "ClusterGraph.hpp"
#include "orderPairs.hpp"
#include "SimilarPairs.hpp"
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

#include "fstream.hpp"
#include "iostream.hpp"
#include "utility.hpp"
#include "vector.hpp"
#include <numeric>
#include <regex>
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
        throw runtime_error("Directory " + directoryName + " already exists.");
    }

    // Create the directory. This guarantees that we start with an empty directory.
    if(!boost::filesystem::create_directory(directoryName)) {
        throw runtime_error("Unable to create directory " + directoryName);
    }

    geneNames.createNew(directoryName + "/" + "GeneNames", parameters.geneCapacity);
    cells.createNew(directoryName + "/" + "Cells");
    cellNames.createNew(directoryName + "/" + "CellNames", parameters.cellCapacity);
    cellMetaData.createNew(directoryName + "/" + "CellMetaData");
    cellMetaDataNames.createNew(directoryName + "/" + "CellMetaDataNames", parameters.cellMetaDataNameCapacity);
    cellMetaDataValues.createNew(directoryName + "/" + "CellMetaDataValues", parameters.cellMetaDataValueCapacity);
    cellMetaDataNamesUsageCount.createNew(directoryName + "/" + "CellMetaDataNamesUsageCount");
    cellExpressionCounts.createNew(directoryName + "/" + "CellExpressionCounts");

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
ExpressionMatrix::ExpressionMatrix(
    const string& directoryName,
    uint64_t geneCapacity,
    uint64_t cellCapacity,
    uint64_t cellMetaDataNameCapacity,
    uint64_t cellMetaDataValueCapacity
    ) :
    ExpressionMatrix(
        directoryName,
        ExpressionMatrixCreationParameters(geneCapacity, cellCapacity, cellMetaDataNameCapacity, cellMetaDataValueCapacity))
{
}



ExpressionMatrixCreationParameters::ExpressionMatrixCreationParameters(
    uint64_t geneCapacity,
    uint64_t cellCapacity,
    uint64_t cellMetaDataNameCapacity,
    uint64_t cellMetaDataValueCapacity
    ) :
    geneCapacity(geneCapacity),
    cellCapacity(cellCapacity),
    cellMetaDataNameCapacity(cellMetaDataNameCapacity),
    cellMetaDataValueCapacity(cellMetaDataValueCapacity)
{
}



// Access a previously created expression matrix stored in the specified directory.
ExpressionMatrix::ExpressionMatrix(const string& directoryName, bool allowReadOnly) :
    directoryName(directoryName)
{
    // Access the binary data with read-write access, so we can add new cells
    // and perform other operations that change the state on disk.
    geneNames.accessExistingReadWrite(directoryName + "/" + "GeneNames", allowReadOnly);
    cells.accessExistingReadWrite(directoryName + "/" + "Cells", allowReadOnly);
    cellNames.accessExistingReadWrite(directoryName + "/" + "CellNames", allowReadOnly);
    cellMetaData.accessExistingReadWrite(directoryName + "/" + "CellMetaData", allowReadOnly);
    cellMetaDataNames.accessExistingReadWrite(directoryName + "/" + "CellMetaDataNames", allowReadOnly);
    cellMetaDataValues.accessExistingReadWrite(directoryName + "/" + "CellMetaDataValues", allowReadOnly);
    cellMetaDataNamesUsageCount.accessExistingReadWrite(directoryName + "/" + "CellMetaDataNamesUsageCount", allowReadOnly);
    cellExpressionCounts.accessExistingReadWrite(directoryName + "/" + "CellExpressionCounts", allowReadOnly);
    cellSets.accessExisting(directoryName, allowReadOnly);
    if(!cellSets.exists("AllCells")) {
        throw runtime_error("Cell set \"AllCells\" is missing.");
    }



    // Access the gene sets.
    using boost::filesystem::directory_iterator;
    std::regex regex(directoryName + "/GeneSet-(.*)-GlobalIds");
    for(auto it = directory_iterator(directoryName); it != directory_iterator(); ++it) {
        const string fileName = it->path().string();
        std::smatch regexMatchResults;
        if(!std::regex_match(fileName, regexMatchResults, regex)) {
            continue;
        }
        CZI_ASSERT(regexMatchResults.size() == 2);
        const string& geneSetName = regexMatchResults[1];
        geneSets[geneSetName].accessExisting(directoryName + "/GeneSet-" + geneSetName, allowReadOnly);
    }
    if(geneSets.find("AllGenes") == geneSets.end()) {
        throw runtime_error("Gene set \"AllGenes\" is missing.");
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
    CZI_ASSERT(geneSets.find("AllGenes") != geneSets.end());

    const StringId stringId = geneNames(geneName);
    if(stringId == geneNames.invalidStringId) {
        const GeneId geneId = GeneId(geneNames[geneName]);
        geneSets["AllGenes"].addGene(geneId);
        geneSets["AllGenes"].forceSorted(); // We guarantee that it remains sorted.
        return true;
    } else {
        return false;   // Was already present.
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
    const vector< pair<string, string> >& metaDataArgument,
    const vector< pair<string, float> >& expressionCounts)
{
#if 0
    cout << "ExpressionMatrix::addCell called." << endl;
    cout << "Cell meta data:" << endl;
    for(const auto& p: metaData) {
        cout << p.first << " " << p.second << endl;
    }
    cout << "Expression counts:" << endl;
    for(const auto& p: expressionCounts) {
        cout << p.first << " " << p.second << endl;
    }
#endif

    // Make a writable copy of the meta data.
    // We will need it to move the cell name to the beginning.
    vector< pair<string, string> > metaData = metaDataArgument;

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
        for(const auto& p : metaData) {
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
        const string& geneName = p.first;
        addGene(geneName);
        const StringId geneId = geneNames(geneName);
        CZI_ASSERT(geneId != geneNames.invalidStringId);
        const float value = p.second;
        if(value < 0.) {
            throw runtime_error("Negative expression count encountered.");
        }
        if(value == 0.) {
            continue;
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
    for(size_t i = 1; i < storedExpressionCounts.size(); i++) {
        if(storedExpressionCounts[i - 1].first == storedExpressionCounts[i].first) {
            const string duplicateGeneName = geneNames[storedExpressionCounts[i].first];
            throw runtime_error("Duplicate expression count for cell " + cellName + " gene " + duplicateGeneName);
        }
    }

    // Add this cell to the AllCells set.
    cellSets.cellSets["AllCells"]->push_back(CellId(cells.size()));

    // Store fixed size information for this cell.
    cells.push_back(cell);


    // Sanity checks.
    CZI_ASSERT(cellNames.size() == cells.size());
    CZI_ASSERT(cellMetaData.size() == cells.size());
    CZI_ASSERT(cellExpressionCounts.size() == cells.size());
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
CellId ExpressionMatrix::addCellFromJson(const string& jsonString)
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
        return addCell(metaData, expressionCounts);

    } catch(...) {

        // If an error occurred, make sure to write out the JSON that caused it.
        cout << "Error processing the following cell JSON:" << endl;
        cout << jsonString << endl;
        throw;
    }

}



// Add cells from data in files with fields separated by commas or by other separators.
// This is a bit messy, but all of this is necessary to allow
// some flexibility in the contents of the input files.
// For example:
// - The first token of the first line of each input file
//   is optional (and ignored is present).
// - Each of the input files can use a different separator.
// - The input files can have Unix or Windows line ends.
// - Expression counts can contain leading and trailing blanks.
// - The ordering of cells is not necessarily the same
//   in the two input files.
// - Some cells can be present in only one input file.
//   Only cells present both in the cell meta data file
//   and in the expression count file are kept.
// See ExpressionMatrix.hpp for usage information.
void ExpressionMatrix::addCells(
    const string& expressionCountsFileName,
    const string& expressionCountsFileSeparators,
    const string& cellMetaDataFileName,
    const string& cellMetaDataFileSeparators
    )
{
    cout << timestamp << "Begin addCells: " << cellCount() <<" cells, "
        << geneCount() << " genes." << endl;

    // Get the number of meta data fields.
    const size_t metaDataCount =
        countTokensInSecondLine(cellMetaDataFileName, cellMetaDataFileSeparators);
    cout << "Cell meta data file " << cellMetaDataFileName
        << " contains " << metaDataCount
        << " meta data fields." << endl;

    // Open the cell meta data file.
    ifstream metaDataFile(cellMetaDataFileName);
    if(!metaDataFile) {
        throw runtime_error("Error opening cell meta data file " + cellMetaDataFileName);
    }

    // Read the header line from the cell meta data file and parse it.
    // This contains meta data names, plus possibly an initial token
    // to be ignored.
    string line;
    getline(metaDataFile, line);
    if(line.empty()) {
        throw runtime_error("Error reading header line from cell meta data file " + cellMetaDataFileName);
    }
    if(line[line.size()-1] == 13) {
        line.resize(line.size()-1); // Remove Windows style line end if necessary.
    }
    vector<string> tokens;
    tokenize(cellMetaDataFileSeparators, line, tokens, false);

    // Get the meta data names.
    vector<string> metaDataNames;
    metaDataNames.push_back("CellName");
    if(tokens.size() == metaDataCount) {
        copy(
            tokens.begin()+1,
            tokens.end(),
            back_inserter(metaDataNames));
    } else if(tokens.size() == metaDataCount-1) {
        copy(
            tokens.begin(),
            tokens.end(),
            back_inserter(metaDataNames));
    } else {
        throw runtime_error("Number of tokens in header line of cell meta data file " +
            cellMetaDataFileName + " is inconsistent with file contents.");
    }
    CZI_ASSERT(metaDataNames.size() == metaDataCount);



    // Read the rest of the cell meta data file.
    // Store the meta data in a map keyed by the cell name.
    vector< vector<string> > metaDataInFile;    // Indexed by CellId within the file
    map<string, CellId> metaDataInFileMap;      // Map from cell name to CellId.
    while(true) {

        // Read a line.
        getline(metaDataFile, line);
        if(!metaDataFile) {
            break;  // We reached the end of the cell meta data file.
        }
        if(line.empty()) {
            throw runtime_error("Empty line in cell meta data file " + cellMetaDataFileName);
        }
        if(line[line.size()-1] == 13) {
            line.resize(line.size()-1); // Remove Windows style line end if necessary.
        }

        // Parse it.
        tokens.clear();
        tokenize(cellMetaDataFileSeparators, line, tokens, false);



        // Check that we have the expected number of tokens.
        if(tokens.size() != metaDataNames.size()) {
            cout << "Unexpected number of items in this line of cell meta data file "
                << cellMetaDataFileName << ":\n";
            cout << line << endl;
            throw runtime_error("Unexpected number of items in line of cell meta data file " + cellMetaDataFileName);
        }


        // Check that we have not already encountered this cell name.
        const string& cellName = tokens[0];
        if(metaDataInFileMap.find(cellName) != metaDataInFileMap.end()) {
            throw runtime_error("Duplicate cell name " + cellName +
                " in cell meta data file " + cellMetaDataFileName);
        }

        // Store the meta data for this cell.
        metaDataInFileMap.insert(make_pair(cellName, metaDataInFile.size()));
        metaDataInFile.push_back(tokens);
        CZI_ASSERT(metaDataInFile.size() == metaDataInFileMap.size());
    }
    CZI_ASSERT(metaDataInFile.size() == metaDataInFileMap.size());
    cout << "Cell meta data file " << cellMetaDataFileName <<
        " contains meta data for " << metaDataInFile.size() << " cells." << endl;



    // Get the number of cells in the expression counts file.
    // We will only keep cells that appear both in the cell meta data file
    // and in the expression counts file.
    const size_t cellCountInExpressionFile =
        countTokensInSecondLine(expressionCountsFileName, expressionCountsFileSeparators) - 1;
    cout << "Cell expression counts file " << expressionCountsFileName <<
        " contains data for " << cellCountInExpressionFile << " cells." << endl;


    // Open the expression counts file.
    ifstream expressionCountsFile(expressionCountsFileName);
    if(!expressionCountsFile) {
        throw runtime_error("Error opening expression counts file " + expressionCountsFileName);
    }


    // Read the names of the cells in the expression counts file.
    getline(expressionCountsFile, line);
    if(line.empty()) {
        throw runtime_error("Error reading header line from expression counts file " + expressionCountsFileName);
    }
    if(line[line.size()-1] == 13) {
        line.resize(line.size()-1); // Remove Windows style line end if necessary.
    }
    tokenize(expressionCountsFileSeparators, line, tokens, true);
    vector<string> cellNamesInExpressionCountsFile;
    if(tokens.size() == cellCountInExpressionFile) {
        cellNamesInExpressionCountsFile = tokens;
    } else if(tokens.size() == cellCountInExpressionFile+1) {
        copy(tokens.begin()+1, tokens.end(), back_inserter(cellNamesInExpressionCountsFile));
    } else {
        throw runtime_error("Number of tokens in header line of expression counts file " +
            expressionCountsFileName + " is inconsistent with file contents.");
    }
    CZI_ASSERT(cellNamesInExpressionCountsFile.size() == cellCountInExpressionFile);



    // We will add the cells in the order used in the meta data file.
    // Create a list of the cells present in both the cell meta data file
    // and the expression counts file. For each store a pair
    // (index in meta data file, index in expression counts file).
    vector< pair<CellId, CellId> > cellsToBeKept;
    for(size_t i=0; i<cellNamesInExpressionCountsFile.size(); i++) {
        const string& cellName = cellNamesInExpressionCountsFile[i];
        const auto it = metaDataInFileMap.find(cellName);
        if(it == metaDataInFileMap.end()) {
            continue;
        }
        const auto j = it->second;
        cellsToBeKept.push_back(make_pair(j, i));
    }
    cout << "The number of cells that appear in both the cell meta data file " <<
        " and the expression counts file is " << cellsToBeKept.size() << "." << endl;
    cout << "This is the number of cells that will be kept." << endl;


    // Gene names encountered in the expression counts file.
    set<string> geneNamesInFile;
    vector<GeneId> geneIdsInFile; // The global gene id, indexed by GeneId in file.

    // Create vectors of pairs (GeneId, count) for each cell to be kept.
    vector< vector< pair<GeneId, float> > > counts(cellsToBeKept.size());



    // Read the rest of the expression counts file.
    while(true) {
        if((geneIdsInFile.size()%1000)==0) {
            cout << timestamp << "Working on gene " << geneIdsInFile.size() << endl;
        }

        // Read a line.
        getline(expressionCountsFile, line);
        if(!expressionCountsFile) {
            break;  // We reached the end of the expression counts file.
        }
        if(line.empty()) {
            throw runtime_error("Empty line in expression counts file " + expressionCountsFileName);
        }
        if(line[line.size()-1] == 13) {
            line.resize(line.size()-1); // Remove Windows style line end if necessary.
        }

        // Parse it.
        tokens.clear();
        tokenize(expressionCountsFileSeparators, line, tokens, true);


        // Check that we have the expected number of tokens.
        if(tokens.size() != cellCountInExpressionFile+1) {
            cout << "Unexpected number of items in this line of expression counts file "
                << expressionCountsFileName << ":\n";
            cout << line << endl;
            throw runtime_error("Unexpected number of items in line of expression counts file " + expressionCountsFileName);
        }

        // Extract the gene name and check that we did not already encounter it in this file.
        const string& geneName = tokens[0];
        if(geneNamesInFile.find(geneName) != geneNamesInFile.end()) {
            throw runtime_error("Duplicate entry for gene " + geneName +
                " in cell expression counts file " + expressionCountsFileName);
        }
        geneNamesInFile.insert(geneName);
        addGene(geneName);
        const GeneId geneId = geneIdFromName(geneName);
        CZI_ASSERT(geneId != invalidGeneId);
        geneIdsInFile.push_back(geneId);



        // Extract the expression counts for this gene,
        // for all the cells we are interested in.
        // For speed in parsing we use strtof instead of boost::lexical_cast,
        // and we check for common representations of zero.
        for(size_t i=0; i<cellsToBeKept.size(); i++) {
            const auto& p = cellsToBeKept[i];
            const CellId cellIdInExpressionCountsFile = p.second;
            const string& countString = tokens[cellIdInExpressionCountsFile+1];
            if(countString=="0" || countString=="0.") {
                continue;
            }
            const char* countStringBegin = countString.c_str();
            char* check = const_cast<char*>(countStringBegin + countString.size());
            const float count = strtof(countStringBegin, &check);
            if(check == countStringBegin) {
                throw runtime_error("Invalid expression count " +
                    (countString.empty() ? "(empty)" : countString) +
                    " encountered at line " +
                    lexical_cast<string>(geneNamesInFile.size()+1) +
                    " of file " + expressionCountsFileName);
            }
            counts[i].push_back(make_pair(geneId, count));
        }
    }



    // Now we can add all the cells to be kept.
    vector< pair<string, string> > metaDataForOneCell;
    vector< pair<string, float> > countsForOneCell;
    for(const string& metaDataName: metaDataNames) {
        metaDataForOneCell.push_back(make_pair(metaDataName, ""));
    }
    for(size_t i=0; i<cellsToBeKept.size(); i++) {
        const auto& p = cellsToBeKept[i];
        const CellId cellIdInMetaDataFile = p.first;
        // const CellId cellIdInExpressionCountsFile = p.second;
        CZI_ASSERT(metaDataInFile[cellIdInMetaDataFile].size() == metaDataNames.size());
        for(size_t j=0; j<metaDataNames.size(); j++) {
            metaDataForOneCell[j].second = metaDataInFile[cellIdInMetaDataFile][j];
        }
        countsForOneCell.clear();
        for(const auto& p: counts[i]) {
            countsForOneCell.push_back(make_pair(geneName(p.first), p.second));
        }
        addCell(metaDataForOneCell, countsForOneCell);
    }

    cout << timestamp << "End addCells: " << cellCount() <<" cells, "
        << geneCount() << " genes." << endl;
}



// Add cells from data in files with fields separated by commas or by other separators.
// See ExpressionMatrix.hpp for usage information.
// This is the old version that needs much more memory.
void ExpressionMatrix::addCellsOld1(
    const string& expressionCountsFileName,
    const string& expressionCountsFileSeparators,
    const string& cellMetaDataFileName,
    const string& cellMetaDataFileSeparators
    )
{
    // Tokenize the cell meta data file and the expression counts file and verify
    // that all lines have the same number of tokens (which cannot be zero), with the possible exception
    // of line one which could have one less token than all other lines.
    vector<vector<string> > cellMetaDataFileLines, expressionCountsFileLines;
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
            CZI_ASSERT(cellMetaDataFileLines[0].size() == cellMetaDataFileLines[1].size() - 1); // This was checked by tokenizeFileAndCheck.
            skip = 0;
        }
        copy(cellMetaDataFileLines[0].begin() + skip, cellMetaDataFileLines[0].end(),
            back_inserter(metaDataFieldNames));
    }



    // Check that there are no duplications in the meta data field names.
    {
        set<string> metaDataFieldNamesSet;
        for(const string& metaDataFieldName : metaDataFieldNames) {
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
            CZI_ASSERT(expressionCountsFileLines[0].size() == expressionCountsFileLines[1].size() - 1); // This was checked by tokenizeFileAndCheck.
            skip = 0;
        }
        copy(expressionCountsFileLines[0].begin() + skip, expressionCountsFileLines[0].end(),
            back_inserter(expressionFileCellNames));
    }
    map<string, CellId> expressionFileCellNamesMap;
    for(size_t i = 0; i < expressionFileCellNames.size(); i++) {
        const string& cellName = expressionFileCellNames[i];
        if(expressionFileCellNamesMap.find(cellName) != expressionFileCellNamesMap.end()) {
            throw runtime_error("Cell " + cellName + " has more than one column in the expression counts file.");
        }
        expressionFileCellNamesMap.insert(make_pair(cellName, i));
    }



    // Summarize the number of cells, genes, and meta data names seen in each file.
    cout << "Cell meta data file " << cellMetaDataFileName << " contains data for ";
    cout << cellMetaDataFileLines.size() - 1 << " cells and ";
    cout << metaDataFieldNames.size() << " meta data names." << endl;
    cout << "Expression counts file " << expressionCountsFileName << " contains data for ";
    cout << expressionFileCellNames.size() << " cells and " << expressionCountsFileLines.size() - 1 << " genes."
        << endl;

    // Add the genes.
    // We want to add them independently of the cells, so they all get added, even the ones
    // for which all cells have zero count.
    CZI_ASSERT(expressionCountsFileLines.size() > 0);   // This was checked by tokenizeFileAndCheck.
    for(size_t i = 1; i < expressionCountsFileLines.size(); i++) {
        const vector<string>& line = expressionCountsFileLines[i];
        CZI_ASSERT(line.size() > 0);    // This was checked by tokenizeFileAndCheck.
        addGene(line.front());
    }



    // Loop over cells in the cell meta data file, but only add the ones that also appear
    // in the expression counts file.
    cout << timestamp << "Storing expression counts and cell meta data." << endl;
    CZI_ASSERT(cellMetaDataFileLines.size() > 1);   // This was checked by tokenizeFileAndCheck.
    CellId addedCellCount = 0;
    for(size_t cellMetaDataFileLine = 1; cellMetaDataFileLine < cellMetaDataFileLines.size(); cellMetaDataFileLine++) {
        if((cellMetaDataFileLine % 1000) == 0) {
            cout << timestamp << "Working on cell meta data file line " << cellMetaDataFileLine + 1 << " of "
                << cellMetaDataFileLines.size() << endl;
        }
        const vector<string>& metaDataLine = cellMetaDataFileLines[cellMetaDataFileLine];
        CZI_ASSERT(metaDataLine.size() > 1);    // This was checked by tokenizeFileAndCheck.
        const string& cellName = metaDataLine.front();

        // See if this cell appears in the expression counts file.
        // If not, skip this cell.
        const auto it = expressionFileCellNamesMap.find(cellName);
        if(it == expressionFileCellNamesMap.end()) {
            continue;   // It's not in the expression counts file. Skip it.
        }

        // Find the column in the expression counts file that contains data for this cell.
        const size_t expressionFileColumn = it->second + 1;

        // Gather the meta data for this cell.
        CZI_ASSERT(metaDataLine.size() == metaDataFieldNames.size() + 1); // This was checked by tokenizeFileAndCheck.
        vector<pair<string, string> > thisCellMetaData;
        thisCellMetaData.push_back(make_pair("CellName", cellName));
        for(size_t i = 0; i < metaDataFieldNames.size(); i++) {
            thisCellMetaData.push_back(make_pair(metaDataFieldNames[i], metaDataLine[i + 1]));
        }

        // Gather the expression counts for this cell.
        vector<pair<string, float> > thisCellExpressionCounts;
        for(size_t i = 1; i < expressionCountsFileLines.size(); i++) {
            const vector<string>& line = expressionCountsFileLines[i];
            CZI_ASSERT(line.size() > 0);    // This was checked by tokenizeFileAndCheck.
            const string& geneName = line.front();
            const string& expressionCountString = line[expressionFileColumn];
            float expressionCount;
            try {
                expressionCount = lexical_cast<float>(expressionCountString);
            } catch (boost::bad_lexical_cast) {
                throw runtime_error("Invalid expression count " + expressionCountString +
                    " for cell " + cellName + " gene " + geneName);
            }
            if(expressionCount != 0.) {
                thisCellExpressionCounts.push_back(make_pair(geneName, expressionCount));
            }
        }

        // Now we can add this cell.
        ++addedCellCount;
        addCell(thisCellMetaData, thisCellExpressionCounts);
    }

    cout << timestamp << "Added " << addedCellCount;
    cout << " cells that appear in both the cell meta data file and the expression counts file." << endl;
    cout << "There are " << cellCount() << " cells and " << geneCount() << " genes." << endl;

}



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



// Return the name of the gene with the given id.
string ExpressionMatrix::geneName(GeneId geneId) const
{
    return geneNames[geneId];
}



// Return a gene id given a string.
// The string can be a gene name or GeneId (a string).
// Returns imnvalidGeneId if the gene was not found.
GeneId ExpressionMatrix::geneIdFromString(const string& s) const
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



// Return a gene id given its name.
// Returns invalidGeneId if the gene was not found.
GeneId ExpressionMatrix::geneIdFromName(const string& geneName) const
{
    return geneNames(geneName);
}



// Return the value of a specified meta data field for a given cell.
// Returns an empty string if the cell does not have the specified meta data field.
string ExpressionMatrix::getCellMetaData(CellId cellId, const string& name) const
{
    // Find the string id of the name.
    // If it does not exist, return an empty string.
    const StringId nameId = cellMetaDataNames(name);
    if(nameId == cellMetaDataNames.invalidStringId) {
        return "";
    }
    return getCellMetaData(cellId, nameId);
}
string ExpressionMatrix::getCellMetaData(CellId cellId, StringId nameId) const
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



// Return a vector containing all of the meta data (Name, Value) pairs
// for a given cell.
vector< pair<string, string> > ExpressionMatrix::getCellMetaData(CellId cellId) const
{
    vector< pair<string, string> > allCellMetaData;
    for(const auto& metaDataPair: cellMetaData[cellId]) {
        const StringId nameId = metaDataPair.first;
        const StringId valueId = metaDataPair.second;
        allCellMetaData.push_back(make_pair(cellMetaDataNames[nameId], cellMetaDataValues[valueId]));
    }
    return allCellMetaData;
}



// Return a vector containing vectors with all of the meta data (Name, Value) pairs
// for a given set of cells.
vector< vector< pair<string, string> > > ExpressionMatrix::getCellMetaData(
    const vector<CellId>& cellIds) const
{
    vector< vector< pair<string, string> > > v(cellIds.size());
    for(size_t i=0; i<cellIds.size(); i++) {
        const CellId cellId = cellIds[i];
        v[i] = getCellMetaData(cellId);
    }
    return v;
}


// Set a meta data (name, value) pair for a given cell.
// If the name already exists for that cell, the value is replaced.
void ExpressionMatrix::setCellMetaData(CellId cellId, const string& name, const string& value)
{
    const StringId nameId = cellMetaDataNames[name];
    const StringId valueId = cellMetaDataValues[value];
    setCellMetaData(cellId, nameId, valueId);
}
void ExpressionMatrix::setCellMetaData(CellId cellId, StringId nameId, const string& value)
{
    const StringId valueId = cellMetaDataValues[value];
    setCellMetaData(cellId, nameId, valueId);
}
void ExpressionMatrix::setCellMetaData(CellId cellId, StringId nameId, StringId valueId)
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



// Get the expression count for a given cell and gene.
// This can be zero.
float ExpressionMatrix::getCellExpressionCount(CellId cellId, GeneId geneId) const
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



// Same as above, specifying a gene name instead of a GeneId.
float ExpressionMatrix::getCellExpressionCount(CellId cellId, const string& geneName) const
{
    // Find the GeneId.
    const GeneId geneId = geneIdFromName(geneName);
    if(geneId == invalidGeneId) {
        throw runtime_error("Gene " + geneName + " does not exist.");
    }

    // Call the function that uses the GeneId.
    return getCellExpressionCount(cellId, geneId);
}



// Get all the non-zero expression counts for a given cell.
vector< pair<GeneId, float> > ExpressionMatrix::getCellExpressionCounts(CellId cellId) const
{
    const auto& countsForThisCell = cellExpressionCounts[cellId];
    vector< pair<GeneId, float> > returnVector;
    returnVector.reserve(countsForThisCell.size());
    for(const auto& p: countsForThisCell) {
        returnVector.push_back(p);
    }
    return returnVector;
}



// Get the expression count for a given gene, for a specified set of cells.
// Each position in the returned vector has the count for
// the cell at the same position in the input vector.
// Some of the returned counts can be zero.
vector<float> ExpressionMatrix::getCellsExpressionCount(const vector<CellId>& cellIds, GeneId geneId) const
{
    vector<float> returnVector;
    returnVector.reserve(cellIds.size());
    for(const CellId cellId: cellIds) {
        returnVector.push_back(getCellExpressionCount(cellId, geneId));
    }
    return returnVector;
}



// Same as above, specifying a gene name instead of a GeneId.
vector<float> ExpressionMatrix::getCellsExpressionCount(const vector<CellId>& cellIds, const string& geneName) const
{
    // Find the GeneId.
    const GeneId geneId = geneIdFromName(geneName);
    if(geneId == invalidGeneId) {
        throw runtime_error("Gene " + geneName + " does not exist.");
    }

    // Call the function that uses the GeneId.
    return getCellsExpressionCount(cellIds, geneId);
}



// Get all the non-zero expression counts for a specified set of cells.
// Each position in the returned vector has the counts for
// the cell at the same position in the input vector.
vector< vector< pair<GeneId, float> > > ExpressionMatrix::getCellsExpressionCounts(const vector<CellId>& cellIds) const
{
    vector< vector< pair<GeneId, float> > > returnVector(cellIds.size());
    for(size_t i=0; i<cellIds.size(); i++) {
        const CellId cellId = cellIds[i];
        vector< pair<GeneId, float> > v = getCellExpressionCounts(cellId);
        v.swap(returnVector[i]);
    }
    return returnVector;
}



// Get all the non-zero expression counts for a specified set of cells,
// but only including a specified set of genes.
// Each position in the returned vector has the counts for
// the cell at the same position in the input vector.
// Note that in each returned pair<GeneId, float>, the geneId is
// a global GeneId.
vector< vector< pair<GeneId, float> > > ExpressionMatrix::getCellsExpressionCountsForGenes(
    const vector<CellId>& cellIds,
    const vector<GeneId>& geneIds) const
{

    // Create a vector of flags that, for each global gene,
    // tells us whether that gene is in the geneIds vector.
    vector<bool> isGeneIncluded(geneCount(), false);
    for(const GeneId& geneId: geneIds) {
        isGeneIncluded[geneId] = true;
    }

    // Create the vector to be returned, with one entry for each of the cell ids.
    vector< vector< pair<GeneId, float> > > returnVector(cellIds.size());



    // Process one cell at a time.
    for(size_t i=0; i<cellIds.size(); i++) {
        const CellId cellId = cellIds[i];

        // Loop over non-zero expression counts for this cell.
        for(const pair<GeneId, float>& p: cellExpressionCounts[cellId]) {

            // Add it, only if this is one of the genes we want.
            const GeneId geneId = p.first;
            if(isGeneIncluded[geneId]) {
                returnVector[i].push_back(p);
            }
        }
    }



    // Done.
    return returnVector;
}



// Compute the average expression vector for a given gene set
// and for a given vector of cells (which is not the same type as a CellSet).
// The last parameter controls the normalization used for the expression counts
// for averaging:
// 0: no normalization (raw read counts).
// 1: L1 normalization (fractional read counts).
// 2: L2 normalization.
void ExpressionMatrix::computeAverageExpression(
    const GeneSet& geneSet,
    const vector<CellId> cellIds,
    vector<double>& averageExpression,
    NormalizationMethod normalizationMethod) const
    {

    // Vector to contain the normalized expression vector for a single cell.
    vector< pair<GeneId, float> > cellExpressionVector;

    // Initialize the average expression to zero.
    averageExpression.resize(geneSet.size());
    fill(averageExpression.begin(), averageExpression.end(), 0.);

    // Accumulate the contribution of all the cells.
    for(const CellId cellId : cellIds) {

        // Compute the normalized expression vector for this cell.
        computeExpressionVector(cellId, geneSet, normalizationMethod, cellExpressionVector);

        // Add all of the expression counts for this cell.
        for(const auto& p : cellExpressionVector) {
            const GeneId localGeneId = p.first;
            const float normalizedCount = p.second;
            averageExpression[localGeneId] += normalizedCount;
        }
    }



    // Divide by the number of cells.
    const double factor = 1. / double(cellIds.size());
    for(double& a : averageExpression) {
        a *= factor;
    }



    // Normalize as requested.
    switch(normalizationMethod) {
    case NormalizationMethod::None:
        break;
    case NormalizationMethod::L1:
        {
        const double factor = 1. / std::accumulate(averageExpression.begin(), averageExpression.end(), 0.);
        for(double& a : averageExpression) {
            a *= factor;
        }
        break;
    }
    case NormalizationMethod::L2:
        {
        double sum = 0.;
        for(const double& a : averageExpression) {
            sum += a * a;
            ;
        }
        const double factor = 1. / sqrt(sum);
        for(double& a : averageExpression) {
            a *= factor;
        }
        break;
    }
    default:
        CZI_ASSERT(0);
    }
}



// Compute the expression vector for a cell and a given GeneSet,
// normalizing it as requested.
// The expression vector contains pairs(local gene id, count).
// The local gene id is an index in the GeneSet.
void ExpressionMatrix::computeExpressionVector(
    CellId cellId,
    const GeneSet& geneSet,
    NormalizationMethod normalizationMethod,
    vector< pair<GeneId, float> >& expressionVector // The computed expression vector.
    ) const
{
    // Copy the expression vector for the cell into the vector passed as an argument.
    expressionVector.clear();
    for(const auto& p: cellExpressionCounts[cellId]) {
        const GeneId globalGeneId = p.first;
        const GeneId localGeneId = geneSet.getLocalGeneId(globalGeneId);
        if(localGeneId != invalidGeneId) {
            expressionVector.push_back(make_pair(localGeneId, p.second));
        }
    }



    // Normalize it as requested.
    float factor;
    double sum = 0.;
    switch(normalizationMethod) {
    case NormalizationMethod::None:
        return;
    case NormalizationMethod::L1:
        for(const auto& p: expressionVector) {
            sum += p.second;
        }
        factor = float(1./sum);
        break;
    case NormalizationMethod::L2:
        for(const auto& p: expressionVector) {
            sum += p.second * p.second;
        }
        factor = float(1./sqrt(sum));
        break;
    default:
        CZI_ASSERT(0);
    }
    for(auto& p: expressionVector) {
        p.second *= factor;
    }
}



// Compute a sorted histogram of a given cell meta data field.
void ExpressionMatrix::histogramMetaData(
    const CellSet& cellSet,
    StringId metaDataNameId,
    vector< pair<string, size_t> >& sortedHistogram) const
{
    // Create the histogram.
    map<string, size_t> histogram;
    for(const CellId cellId: cellSet) {
        const string metaDataValue = getCellMetaData(cellId, metaDataNameId);
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
// expression counts. This takes into account all genes.
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

#if 0
    cout << "ExpressionMatrix::computeCellSimilarity" << endl;
    cout << "n " << n << endl;
    cout << "scalarProduct " << scalarProduct << endl;
    cout << "cell0.sum1 " << cell0.sum1 << endl;
    cout << "cell1.sum1 " << cell1.sum1 << endl;
    cout << "cell0.sum2 " << cell0.sum2 << endl;
    cout << "cell1.sum2 " << cell1.sum2 << endl;
    cout << "numerator " << numerator << endl;
    cout << "denominator " << denominator << endl;
#endif

    return numerator / denominator;
}


// Compute the similarity between two cells given their CellId.
// The similarity is the correlation coefficient of their
// expression counts. This takes into account
// only the genes in the specified gene set.
// This is done without creating an ExpressionMatrixSubset.
double ExpressionMatrix::computeCellSimilarity(
    const string& geneSetName,
    CellId cellId0,
    CellId cellId1) const
{
    const auto it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    const GeneSet& geneSet = it->second;
    return computeCellSimilarity(geneSet, cellId0, cellId1);
}
double ExpressionMatrix::computeCellSimilarity(
    const GeneSet& geneSet,
    CellId cellId0,
    CellId cellId1) const
{
    // Compute the scalar product of the expression counts for the two cells,
    // taking into account only genes in this gene set.
    typedef pair<GeneId, float>const* Iterator;
    const Iterator begin0 = cellExpressionCounts.begin(cellId0);
    const Iterator end0 = cellExpressionCounts.end(cellId0);
    const Iterator begin1 = cellExpressionCounts.begin(cellId1);
    const Iterator end1 = cellExpressionCounts.end(cellId1);
    Iterator it0 = begin0;
    Iterator it1 = begin1;
    double sum01 = 0.;
    while((it0 != end0) && (it1 != end1)) {
        const GeneId geneId0 = it0->first;
        const GeneId geneId1 = it1->first;

        if(geneId0 < geneId1) {
            ++it0;
        } else if(geneId1 < geneId0) {
            ++it1;
        } else {
            // Here, geneId0==geneId1.
            if(geneSet.contains(geneId0)) {
                sum01 += it0->second * it1->second;
            }
            ++it0;
            ++it1;
        }
    }



    // Compute the other sums we need to compute the correlation coefficient.
    double sum0 = 0.;
    double sum00 = 0.;
    for(it0=begin0; it0!=end0; ++it0) {
        const GeneId geneId0 = it0->first;
        if(geneSet.contains(geneId0)) {
            const auto& count = it0->second;
            sum0 += count;
            sum00 += count*count;
        }
    }
    double sum1 = 0.;
    double sum11 = 0.;
    for(it1=begin1; it1!=end1; ++it1) {
        const GeneId geneId1 = it1->first;
        if(geneSet.contains(geneId1)) {
            const auto& count = it1->second;
            sum1 += count;
            sum11 += count*count;
        }
    }



    // Compute the correlation coefficient.
    // See, for example, https://en.wikipedia.org/wiki/Correlation_and_dependence
    const double n = geneSet.size();
    const double numerator = n*sum01 - sum0*sum1;
    const double denominator = sqrt(
            (n*sum00 - sum0*sum0) *
            (n*sum11 - sum1*sum1)
            );
    return numerator / denominator;

}





// Get the names of all currently defined cell sets.
vector<string> ExpressionMatrix::getCellSetNames() const
{
    vector<string> names;
    for(const auto& p: cellSets.cellSets) {
        names.push_back(p.first);
    }
    return names;
}



// Create a new cell set that contains cells for which
// the value of a specified meta data field is identical to a given string
// or matches a given regular expression.
// Return true if successful, false if a cell set with
// the specified name already exists.
bool ExpressionMatrix::createCellSetUsingMetaData(
    const string& cellSetName,          // The name of the cell set to be created.
    const string& metaDataFieldName,    // The name of the meta data field to be used.
    const string& matchString,          // The string or regular expression that must be matched.
    bool useRegex                       // true=match as regular expression, false=match as string.
    )
{
    // See if a cell set with this name already exists.
    if(cellSets.exists(cellSetName)) {
        cout << "Cell set " << cellSetName << " already exists." << endl;
        return false;
    }

    // If regular expression matching was requested, create the regular expression we are going to match.
    std::regex regex;
    if(useRegex) {
        regex = matchString;
    }



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

            // Figure out if this cell should be included in the new cell set.
            const auto metaDataValue = cellMetaDataValues(p.second);
            bool includeThisCell;
            if(useRegex) {
                includeThisCell = std::regex_match(metaDataValue.begin(), metaDataValue.end(), regex);
            } else {
                includeThisCell =
                    (metaDataValue.size() == matchString.size()) &&
                    equal(matchString.begin(), matchString.end(), metaDataValue.begin());
            }

            // If the meta data value matches the given string or regular expression,
            // add the cell to the set.
            if(includeThisCell) {
                cellSet.push_back(cellId);
            }

            // We don't need to check the remaining meta data fields,
            // since we found the one we were looking for.
            break;
        }
    }



    // Store this cell set.
    cellSets.addCellSet(cellSetName, cellSet);
    // cout << "New cell set " << cellSetName << " contains " << cellSet.size() << " cells." << endl;

    return true;
}

// Create a new cell set using a vector of cellIds.
bool ExpressionMatrix::createCellSet(const string& cellSetName, vector<CellId>& cellIds)
{
    if(cellSets.exists(cellSetName)) {
        cout << "Cell set " << cellSetName << " already exists." << endl;
        return false;
    }
    cellSets.addCellSet(cellSetName, cellIds);
    return true;
};



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
    const CellSet& inputSet0 = *(it0->second);
    const auto it1 = cellSets.cellSets.find(inputSetName1);
    if(it1 == cellSets.cellSets.end()) {
        cout << "Cell set " << inputSetName1 << " does not exists." << endl;
        return false;
    }
    const CellSet& inputSet1 = *(it1->second);



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
    const CellSet& inputCellSet = *(it->second);

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
    for(const CellId cellId : inputCellSet) {
        if(uniformGenerator() < probability) {
            outputCellSet.push_back(cellId);
        }
    }

    // Store the new cell set.
    cellSets.addCellSet(outputCellSetName, outputCellSet);

    return true;
}



// Get the names of all currently defined cell similarity graphs.
vector<string> ExpressionMatrix::getCellGraphNames() const
{
    vector<string> v;
    for(const auto& p: cellGraphs) {
        v.push_back(p.first);
    }
    return v;
}



// Create a new graph.
// Graphs are not persistent (they are stored in memory only).
void ExpressionMatrix::createCellGraph(
    const string& graphName,            // The name of the graph to be created. This is used as a key in the graph map.
    const string& cellSetName,          // The cell set to be used.
    const string& similarPairsName,     // The name of the SimilarPairs object to be used to create the graph.
    double similarityThreshold,         // The minimum similarity to create an edge.
    size_t maxConnectivity              // The maximum number of neighbors (k of the k-NN graph).
 )
{
    // A graph with this name should not already exist.
    if(cellGraphs.find(graphName) != cellGraphs.end()) {
        throw runtime_error("Graph " + graphName + " already exists.");
    }

    // Locate the cell set.
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exists.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);

    // Create the graph.
    typedef boost::shared_ptr<CellGraph> GraphSharedPointer;
    const GraphSharedPointer graph = GraphSharedPointer(new CellGraph(
        cellSet,
        directoryName + "/SimilarPairs-" + similarPairsName,
        similarityThreshold,
        maxConnectivity
        ));

    // Create the GraphInformation object that will be stored with the graph.
    CellGraphInformation graphInformation;
    graphInformation.cellSetName = cellSetName;
    graphInformation.similarPairsName = similarPairsName;
    graphInformation.similarityThreshold = similarityThreshold;
    graphInformation.maxConnectivity = maxConnectivity;

    // Remove isolated vertices.
    graphInformation.isolatedVertexCount = graph->removeIsolatedVertices();
    graphInformation.vertexCount = num_vertices(*graph);
    graphInformation.edgeCount = num_edges(*graph);

    // Store it.
    cellGraphs.insert(make_pair(graphName, make_pair(graphInformation, graph)));

}



// Compute the layout (vertex positions) for the graph with a given name.
void ExpressionMatrix::computeCellGraphLayout(const string& graphName)
{
    // Locate the graph.
    const auto it = cellGraphs.find(graphName);
    if(it == cellGraphs.end()) {
        throw runtime_error("Graph " + graphName + " does not exist.");
    }
    CellGraph& cellGraph = *(it->second.second);

    if(!cellGraph.layoutWasComputed) {
        cellGraph.computeLayout();
        cellGraph.layoutWasComputed = true;
    }

}



// Return vertex information for the graph with a given name.
vector<CellGraphVertexInfo> ExpressionMatrix::getCellGraphVertices(const string& graphName) const
{
    // Locate the graph.
    const auto it = cellGraphs.find(graphName);
    if(it == cellGraphs.end()) {
        throw runtime_error("Graph " + graphName + " does not exist.");
    }
    const CellGraph& cellGraph = *(it->second.second);

    if(!cellGraph.layoutWasComputed) {
        throw runtime_error("Layout for graph " + graphName + " is not available.");
    }

    // Fill the return vector by looping over all vertices.
    vector<CellGraphVertexInfo> vertexInfos;
    BGL_FORALL_VERTICES(v, cellGraph, CellGraph) {
        vertexInfos.push_back(cellGraph[v]);
    }
    return vertexInfos;
}



// Return the cell ids of the two vertices corresponding to
// each of the edges of the cell graph with given name.
vector< pair<CellId, CellId> > ExpressionMatrix::getCellGraphEdges(const string& graphName) const
{
    // Locate the graph.
    const auto it = cellGraphs.find(graphName);
    if(it == cellGraphs.end()) {
        throw runtime_error("Graph " + graphName + " does not exist.");
    }
    const CellGraph& cellGraph = *(it->second.second);

    // Loop over graph edges.
    vector< pair<CellId, CellId> > v;
    BGL_FORALL_EDGES(e, cellGraph, CellGraph) {
        const CellGraph::vertex_descriptor v0 = source(e, cellGraph);
        const CellGraph::vertex_descriptor v1 = target(e, cellGraph);
        const CellGraphVertex& vertex0 = cellGraph[v0];
        const CellGraphVertex& vertex1 = cellGraph[v1];
        v.push_back(make_pair(vertex0.cellId, vertex1.cellId));
    }
    return v;
}



// Store the cluster ids in a graph in a meta data field.
void ExpressionMatrix::storeClusterId(
    const string& metaDataName,
    const CellGraph& graph)
{
    // Find the string id corresponding to the specified meta data name.
    // This adds it to the table if not already present.
    const StringId metaDataNameStringId = cellMetaDataNames[metaDataName];

    // Loop over all vertices in the graph.
    BGL_FORALL_VERTICES(v, graph, CellGraph) {
        const CellGraphVertex& vertex = graph[v];

        // Extract the cell id and the cluster id.
        const CellId cellId = vertex.cellId;
        const uint32_t clusterId = vertex.clusterId;

        // Store the cluster id as cell meta data.
        // If the name already exists for this cell, the value is replaced.
        setCellMetaData(cellId, metaDataNameStringId, lexical_cast<string>(clusterId));
    }
}



// Compute gene information content in bits for a given gene set and cell set,
// using the specified normalization method.
// We do it one gene at a time to avoid the need for an amount of
// memory proportional to the product of the number of cells
// times the number of genes.
void ExpressionMatrix::computeGeneInformationContent(
    const GeneSet& geneSet,
    const CellSet& cellSet,
    NormalizationMethod normalizationMethod,
    vector<float>& geneInformationContent) const
    {
    geneInformationContent.reserve(geneSet.size());
    geneInformationContent.clear();
    for(const GeneId geneId : geneSet) {
        geneInformationContent.push_back(computeGeneInformationContent(geneId, cellSet, normalizationMethod));
    }

}



float ExpressionMatrix::computeGeneInformationContent(
    GeneId geneId,
    const CellSet& cellSet,
    NormalizationMethod normalizationMethod) const
{

    // Create a vector of expression counts for this gene and for all cells in the cell set,
    // using the requested normalization.
    // Note that we use the normalization defined using all genes.
    vector<float> count;
    count.reserve(cellSet.size());
    for(const CellId cellId : cellSet) {
        const Cell& cell = cells[cellId];
        float c = getCellExpressionCount(cellId, geneId);
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
    for(const float c : count) {
        sum += double(c);
    }

    // Compute the information content.
    double informationContent = log(double(cellSet.size())); // Equally distributed.
    const double inverseSum = 1. / sum; // No problem with division by zero - never used if sum is zero
    for(const float c : count) {
        if(c > 0.) {
            const double p = c * inverseSum;
            informationContent += p * log(p);
        }
    }


    // Convert to bits.
    informationContent /= log(2.);

    return float(informationContent);
}



void ExpressionMatrix::createGeneSetUsingInformationContent(
    const string& existingGeneSetName,
    const string& cellSetName,
    NormalizationMethod normalizationMethod,
    double geneInformationContentThreshold,
    const string& newGeneSetName)
{
    // Locate the existing gene set.
    const auto itExistingGeneSet = geneSets.find(existingGeneSetName);
    if(itExistingGeneSet == geneSets.end()) {
        throw runtime_error("Gene set " + existingGeneSetName + " does not exist.");
    }
    const GeneSet& existingGeneSet = itExistingGeneSet->second;

    // Locate the cell set.
    const auto itCellSet = cellSets.cellSets.find(cellSetName);
    if(itCellSet == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const CellSet& cellSet = *(itCellSet->second);

    // Verify that the new gene set does not already exist.
    if(geneSets.find(newGeneSetName) != geneSets.end()) {
        throw runtime_error("Gene set " + newGeneSetName + " already exists.");
    }

    // Create the new gene set.
    GeneSet& newGeneSet = geneSets[newGeneSetName];
    newGeneSet.createNew(directoryName + "/GeneSet-" + newGeneSetName);
    createGeneSetUsingInformationContent(
        existingGeneSet,
        cellSet,
        normalizationMethod,
        geneInformationContentThreshold,
        newGeneSet);
}



void ExpressionMatrix::createGeneSetUsingInformationContent(
    const GeneSet& existingGeneSet,
    const CellSet& cellSet,
    NormalizationMethod normalizationMethod,
    double geneInformationContentThreshold,
    GeneSet& newGeneSet) const
{
    // Check that we are starting with an empty set.
    CZI_ASSERT(newGeneSet.size() == 0);

    // Compute gene information content using the requested normalization method.
    vector<float> informationContent;
    computeGeneInformationContent(existingGeneSet, cellSet, normalizationMethod, informationContent);

    // Add to the new gene set the genes that have sufficient information.
    for(GeneId localGeneId=0; localGeneId!=existingGeneSet.size(); localGeneId++) {
        if(informationContent[localGeneId] > geneInformationContentThreshold) {
            const GeneId globalGeneId = existingGeneSet.getGlobalGeneId(localGeneId);
            newGeneSet.addGene(globalGeneId);
        }
    }
}



// Create a new named ClusterGraph by running clustering on an existing CellGraph.
void ExpressionMatrix::createClusterGraph(
    const string& cellGraphName,            // The name of the cell graph to do clustering on.
    const string& clusterGraphName,         // The name of the ClusterGraph to be created.
    size_t stableIterationCount,            // Stop after this many iterations without changes.
    size_t maxIterationCount,               // Stop after this many iterations no matter what.
    size_t seed,                            // To initialize label propagation algorithm.
    size_t minClusterSize,                  // Minimum number of cells for a cluster to be retained.
    size_t maxConnectivity,
    double similarityThreshold              // For edges of the cluster graph.
 )
{
    ClusterGraphCreationParameters parameters(
        stableIterationCount, maxIterationCount, seed, minClusterSize, maxConnectivity, similarityThreshold);
    createClusterGraph(cellGraphName, parameters, clusterGraphName);
}
void ExpressionMatrix::createClusterGraph(
    const string& cellGraphName,
    const ClusterGraphCreationParameters& clusterGraphCreationParameters,
    const string& clusterGraphName)
{
    createClusterGraph(cout, cellGraphName, clusterGraphCreationParameters, clusterGraphName);
}
void ExpressionMatrix::createClusterGraph(
    ostream& out,
    const string& cellGraphName,
    const ClusterGraphCreationParameters& clusterGraphCreationParameters,
    const string& clusterGraphName)
{
    // Locate the cell graph.
    const auto it = cellGraphs.find(cellGraphName);
    if(it == cellGraphs.end()) {
        throw runtime_error("Cell graph " + cellGraphName + " does not exist.");
        return;
    }
    const CellGraphInformation& cellGraphInformation = it->second.first;
    const string& similarPairsName = cellGraphInformation.similarPairsName;
    const SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, true);
    const GeneSet& geneSet = similarPairs.getGeneSet();
    CellGraph& cellGraph = *(it->second.second);



    // Check that a ClusterGraph with this name does not already exists.
    if(clusterGraphs.find(clusterGraphName) != clusterGraphs.end()) {
        throw runtime_error("Cluster graph " + clusterGraphName + " already exists.");
        return;
    }



    // Do the clustering on this cell graph, using the specified parameters.
    cellGraph.labelPropagationClustering(
        out,
        clusterGraphCreationParameters.seed,
        clusterGraphCreationParameters.stableIterationCount,
        clusterGraphCreationParameters.maxIterationCount);



    // Create the ClusterGraph.
    const boost::shared_ptr<ClusterGraph> clusterGraphPointer =
        boost::shared_ptr<ClusterGraph>(new ClusterGraph(cellGraph, geneSet));
    clusterGraphs.insert(make_pair(clusterGraphName, clusterGraphPointer));
    ClusterGraph& clusterGraph = *clusterGraphPointer;

    // Remove the vertices that correspond to small clusters.
    clusterGraph.removeSmallVertices(clusterGraphCreationParameters.minClusterSize);

    // Compute the average expression for each cluster - that is, for each vertex
    // of the cluster graph.
    BGL_FORALL_VERTICES(v, clusterGraph, ClusterGraph) {
        ClusterGraphVertex& vertex = clusterGraph[v];
        const NormalizationMethod normalizationMethod = NormalizationMethod::L2;    // Use L2 normalization. We might need to make this configurable.
        computeAverageExpression(
            geneSet,
            vertex.cells,
            vertex.averageGeneExpression,
            normalizationMethod);
    }

    // Store in each edge the similarity of the two clusters, computed using the clusters
    // average expression stored in each vertex.
    clusterGraph.computeSimilarities();

    // Remove edges with low similarity.
    clusterGraph.removeWeakEdges(clusterGraphCreationParameters.similarityThreshold);

    // Make it a k-nn graph.
    clusterGraph.makeKnn(clusterGraphCreationParameters.maxConnectivity);

    out << "Cluster graph " << clusterGraphName << " has " << num_vertices(clusterGraph);
    out << " vertices and " << num_edges(clusterGraph) << " edges." << endl;
}



// Compute svg and pdf layout with labels for a named cluster graph.
void ExpressionMatrix::computeClusterGraphLayout(
    const string& clusterGraphName,
    size_t timeoutSeconds,
    bool withLabels)
{
    // Locate the cluster graph.
    const auto it = clusterGraphs.find(clusterGraphName);
    if(it == clusterGraphs.end()) {
        throw runtime_error("Cluster graph " + clusterGraphName + " does not exist.");
    }
    ClusterGraph& clusterGraph = *(it->second);

    // Compute the layout.
    clusterGraph.computeLayout(timeoutSeconds, clusterGraphName, geneNames, withLabels);

}



// Get a vector of cluster ids for the vertices of a named cluster graph.
vector<uint32_t> ExpressionMatrix::getClusterGraphVertices(const string& clusterGraphName) const
{
    // Locate the cluster graph.
    const auto it = clusterGraphs.find(clusterGraphName);
    if(it == clusterGraphs.end()) {
        throw runtime_error("Cluster graph " + clusterGraphName + " does not exist.");
    }
    ClusterGraph& clusterGraph = *(it->second);

    // Gather the cluster ids of the vertices.
    vector<uint32_t> clusterIds;
    CZI_ASSERT(clusterGraph.vertexMap.size() == num_vertices(clusterGraph));
    for(const auto& p: clusterGraph.vertexMap) {
        clusterIds.push_back(p.first);
    }

    return clusterIds;
}



// Get the global GeneId's of the genes used by a ClusterGraph.
vector<GeneId> ExpressionMatrix::getClusterGraphGenes(const string& clusterGraphName) const
{
    // Locate the cluster graph.
    const auto it = clusterGraphs.find(clusterGraphName);
    if(it == clusterGraphs.end()) {
        throw runtime_error("Cluster graph " + clusterGraphName + " does not exist.");
    }
    ClusterGraph& clusterGraph = *(it->second);

    // The genes are stored in the ClusterGraph.
    return clusterGraph.geneSet;

}



// Get a vector of the cell ids in a given cluster.
vector<CellId> ExpressionMatrix::getClusterCells(
    const string& clusterGraphName,
    uint32_t clusterId) const
{
    // Locate the cluster graph.
    const auto it = clusterGraphs.find(clusterGraphName);
    if(it == clusterGraphs.end()) {
        throw runtime_error("Cluster graph " + clusterGraphName + " does not exist.");
    }
    ClusterGraph& clusterGraph = *(it->second);

    // Find the vertex corresponding to the requested cluster id.
    const auto jt = clusterGraph.vertexMap.find(clusterId);
    if(jt == clusterGraph.vertexMap.end()) {
        throw runtime_error("Cluster " + lexical_cast<string>(clusterId) +
            " of cluster graph " + clusterGraphName + " does not exist.");
    }
    const ClusterGraph::vertex_descriptor v = jt->second;
    const ClusterGraphVertex& vertex = clusterGraph[v];
    CZI_ASSERT(vertex.clusterId == clusterId);

    // The cells are stored in the vertex.
    return vertex.cells;
}



// Get the average expression vector(L2-normalized) in a cluster.
// The entries in the returned vector correspond on-on-one
// to the gene ids returned by getClusterGraphGenes.
vector<double> ExpressionMatrix::getClusterAverageExpression(
    const string& clusterGraphName,
    uint32_t clusterId) const
{
    // Locate the cluster graph.
    const auto it = clusterGraphs.find(clusterGraphName);
    if(it == clusterGraphs.end()) {
        throw runtime_error("Cluster graph " + clusterGraphName + " does not exist.");
    }
    ClusterGraph& clusterGraph = *(it->second);

    // Find the vertex corresponding to the requested cluster id.
    const auto jt = clusterGraph.vertexMap.find(clusterId);
    if(jt == clusterGraph.vertexMap.end()) {
        throw runtime_error("Cluster " + lexical_cast<string>(clusterId) +
            " of cluster graph " + clusterGraphName + " does not exist.");
    }
    const ClusterGraph::vertex_descriptor v = jt->second;
    const ClusterGraphVertex& vertex = clusterGraph[v];
    CZI_ASSERT(vertex.clusterId == clusterId);

    // The average expression is stored in the vertex.
    return vertex.averageGeneExpression;

}



// Create meta data from the cluster ids stored in a ClusterGraph.
void ExpressionMatrix::createMetaDataFromClusterGraph(
    const string& clusterGraphName,
    const string& metaDataName)
{
    // Locate the cluster graph.
    const auto it = clusterGraphs.find(clusterGraphName);
    if(it == clusterGraphs.end()) {
        throw runtime_error("Cluster graph " + clusterGraphName + " does not exist.");
    }
    ClusterGraph& clusterGraph = *(it->second);

    // Loop over all vertices of the cluster graph.
    // Each vertex corresponds to a cluster.
    BGL_FORALL_VERTICES(v, clusterGraph, ClusterGraph) {
        const ClusterGraphVertex& vertex = clusterGraph[v];
        const string clusterIdString = lexical_cast<string>(vertex.clusterId);
        for(const CellId cellId: vertex.cells) {
            setCellMetaData(cellId, metaDataName, clusterIdString);
        }
    }

    // Loop over unclustered cells.
    for(const CellId cellId: clusterGraph.unclusteredCells) {
        setCellMetaData(cellId, metaDataName, "Unclustered-" + lexical_cast<string>(cellId));
    }
}
