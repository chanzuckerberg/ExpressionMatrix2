#include "ExpressionMatrix.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;


// Return a reference to the cell set with a given name,
// throwing an exception if it does not exist.
CellSet& ExpressionMatrix::cellSet(const string& cellSetName)
{
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    } else {
        return *(it->second);
    }
}
const CellSet& ExpressionMatrix::cellSet(const string& cellSetName) const
{
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    } else {
        return *(it->second);
    }
}

// Check that a cell set with the given name does not exist,
// and throw an exception if it does exist.
void ExpressionMatrix::checkCellSetDoesNotExist(const string& cellSetName) const
{
    if(cellSets.cellSets.find(cellSetName) != cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " already exists.");
    }

}



// Create a new cell set consisting of cells for which a given meta data field
// is numeric and is greater than, less than, or between specified values.
void ExpressionMatrix::createCellSetUsingNumericMetaDataGreaterThan(
    const string& cellSetName,          // The name of the cell set to be created.
    const string& metaDataFieldName,
    double lowerBound)
{
    createCellSetUsingNumericMetaData(cellSetName, metaDataFieldName, true, lowerBound, false, 0.);
}
void ExpressionMatrix::createCellSetUsingNumericMetaDataLessThan(
    const string& cellSetName,          // The name of the cell set to be created.
    const string& metaDataFieldName,
    double upperBound)
{
    createCellSetUsingNumericMetaData(cellSetName, metaDataFieldName, false, 0., true, upperBound);
}
void ExpressionMatrix::createCellSetUsingNumericMetaDataBetween(
    const string& cellSetName,          // The name of the cell set to be created.
    const string& metaDataFieldName,
    double lowerBound,
    double upperBound)
{
    createCellSetUsingNumericMetaData(cellSetName, metaDataFieldName, true, lowerBound, true, upperBound);
}
void ExpressionMatrix::createCellSetUsingNumericMetaData(
    const string& cellSetName,          // The name of the cell set to be created.
    const string& metaDataFieldName,
    bool useLowerBound, double lowerBound,
    bool useUpperBound, double upperBound)
{
    // Check that a cell set with this name does not exist.
    checkCellSetDoesNotExist(cellSetName);

    // Find the cells that belong to the new cell set.
    vector<CellId> cellSet;
    for(CellId cellId=0; cellId<cells.size(); cellId++) {

        // Loop over the meta data fields for this cell.
        const auto metaData = cellMetaData[cellId];
        for(const pair<StringId, StringId>& p: metaData) {

            // If the meta data name is not the specified one, skip.
            if(cellMetaDataNames.equal(p.first, metaDataFieldName)) {
                // We found the meta data field we are looking for. Get its value.
                const string metaDataValue = cellMetaDataValues[p.second];

                // If the value is numeric, check it against the specified bounds.
                try {
                    const double value = lexical_cast<double>(metaDataValue);
                    const bool lowerBoundViolated = useLowerBound && value<lowerBound;
                    const bool upperBoundViolated = useUpperBound && value>upperBound;
                    if(!(lowerBoundViolated || upperBoundViolated)) {
                        cellSet.push_back(cellId);
                    }
                } catch(bad_lexical_cast) {
                    // It was not a number. Don't add this cell to the cell set.
                }

                // We don't need to check the remaining meta data fields for this cell.
                break;
            }

        }
    }



    // Store this cell set.
    cellSets.addCellSet(cellSetName, cellSet);
}



// Remove an existing cell set.
void ExpressionMatrix::removeCellSet(const string& cellSetName)
{
    cellSets.removeCellSet(cellSetName);
}
