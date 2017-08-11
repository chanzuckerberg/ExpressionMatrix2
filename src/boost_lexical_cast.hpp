#ifndef CZI_EXPRESSION_MATRIX2_BOOST_LEXICAL_CAST_HPP
#define CZI_EXPRESSION_MATRIX2_BOOST_LEXICAL_CAST_HPP

// Add to the ExpressionMatrix2 namespace some names defined in namespace boost.

#include <boost/lexical_cast.hpp>

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
    using boost::lexical_cast;
    using boost::bad_lexical_cast;
    }
}

#endif
