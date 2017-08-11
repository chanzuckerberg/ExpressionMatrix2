#ifndef CZI_EXPRESSION_MATRIX2_TOUCH_MEMORY_HPP
#define CZI_EXPRESSION_MATRIX2_TOUCH_MEMORY_HPP

#include "cstddef.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        // Touch a range of memory in order to cause the
        // supporting pages of virtual memory to be loaded in real memory.
        // The return value can be ignored.
        size_t touchMemory(const void* begin, const void* end, size_t pageSize=4096);
    }
}

#endif
