#include "CellSets.hpp"
#include "deduplicate.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>



// Creare a new CellSets in the specified directory.
// We just remove any preexisting CellSets files.
void CellSets::createNew(const string& directoryNameArgument)
{
    // Store the directory name.
    // It will be used if new cell sets are created.
    directoryName = directoryNameArgument;

    // Start with no cell sets.
    cellSets.clear();

    // Loop over all files in the directory.
    boost::regex regex(directoryName + "/CellSet-.*");
    using boost::filesystem::directory_iterator;
    for(auto it=directory_iterator(directoryName); it!=directory_iterator(); ++it) {
        const string fileName = it->path().string();
        if(boost::regex_match(fileName, regex)) {
            boost::filesystem::remove(fileName);
        }
    }
}

// Access existing CellSets in the specified directory.
void CellSets::accessExisting(const string& directoryNameArgument)
{
    // Some boost functionality we usein this function.
    using boost::filesystem::directory_iterator;
    using boost::filesystem::directory_iterator;

    // Store the directory name.
    // It will be used if new cell sets are created.
    directoryName = directoryNameArgument;

    // Start with no cell sets.
    cellSets.clear();

    // Loop over all files in the directory.
    boost::regex regex(directoryName + "/CellSet-.*");
    for(auto it=directory_iterator(directoryName); it!=directory_iterator(); ++it) {
        const string fileName = it->path().string();
        if(!boost::regex_match(fileName, regex)) {
            continue;
        }

        // We found a file containing a CellSet. Access it.
        boost::shared_ptr<CellSet> mappedCellSet = boost::shared_ptr<CellSet>(new CellSet);
        mappedCellSet->accessExistingReadWrite(fileName);

        // Store it in our table of known cell sets.
        const string cellSetName = fileName.substr((directoryName + "/CellSet-").size());
        cellSets.insert(make_pair(cellSetName, mappedCellSet));
    }

}



// Add a new cell set.
void CellSets::addCellSet(
    const string& cellSetName,
    vector<CellId>& cellSet)
{
    // Deduplicate the given cell set.
    ExpressionMatrix2::deduplicate(cellSet);

    // Create the new cell set in mapped memory.
    boost::shared_ptr<CellSet> mappedCellSet = boost::shared_ptr<CellSet>(new CellSet);
    mappedCellSet->createNew(directoryName + "/CellSet-" + cellSetName, cellSet.size());

    // Copy the data from the given cell set.
    copy(cellSet.begin(), cellSet.end(), mappedCellSet->begin());

    // Store it in our table of known cell sets.
    cellSets.insert(make_pair(cellSetName, mappedCellSet));

}
