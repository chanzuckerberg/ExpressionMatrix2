#ifndef CZI_EXPRESSION_MATRIX2_CELL_SETS_HPP
#define CZI_EXPRESSION_MATRIX2_CELL_SETS_HPP

#include "Ids.hpp"
#include "MemoryMappedVector.hpp"

#include <boost/shared_ptr.hpp>

#include "map.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class CellSets;
        using CellSet = MemoryMapped::Vector<CellId>;
    }
}


// Class to maintain sets of cells.
// Each set of cells is a sorted MemoryMapped::Vector<CellId>
// without duplications.

class ChanZuckerberg::ExpressionMatrix2::CellSets {
public:

    // Create a new CellSets in the specified directory.
    void createNew(const string& directoryName);

    // Access existing CellSets in the specified directory.
    void accessExisting(const string& directoryName);

    // Add a new cell set.
    void addCellSet(
        const string& cellSetName,
        vector<CellId>& cellSet);

    // Add a cell to an existing cell set.
    // If the last argument is true, just add the given cell id at the end
    // without checking for duplications and without resorting.
    // The caller is responsible for making sure that the set stays sorted
    // and without duplications.
    void addCell(
        const string& cellSetName,
        CellId,
        bool addAtEnd = false
        );

    // Resort and deduplicate a given set.
    void deduplicate(const string& cellSetName);

    // Return true if a cell set with a given name exists.
    bool exists(const string& cellSetName) const
    {
        return cellSets.find(cellSetName) != cellSets.end();
    }

    // The currently defined cell sets.
    map<string, boost::shared_ptr<CellSet> > cellSets;

    // The name of the directory where files for these CellSets are.
    // The names are patterned after directoryName/CellSet-Name.
    string directoryName;
};


#endif
