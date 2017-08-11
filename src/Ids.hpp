#ifndef CZI_EXPRESSION_MATRIX2_IDS_HPP
#define CZI_EXPRESSION_MATRIX2_IDS_HPP

#include "cstdint.hpp"
#include <limits>

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        // Integers typed use as gene and cell ids.
        // These limit the number of genes and cells that can be used.
        using GeneId = uint32_t;
        using CellId = uint32_t;
        static const GeneId invalidGeneId = std::numeric_limits<GeneId>::max();
        static const CellId invalidCellId = std::numeric_limits<CellId>::max();

        // Type used to identify strings used for cell meta data and for gene names.
        using StringId = uint32_t;
    }
};

#endif
