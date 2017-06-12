#ifndef CZI_EXPRESSION_MATRIX2_TIME_STAMP_HPP
#define CZI_EXPRESSION_MATRIX2_TIME_STAMP_HPP


#include "boost/date_time/posix_time/posix_time.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        inline ostream& timestamp(ostream& s)
        {
            s << boost::posix_time::microsec_clock::local_time() << " ";
            return s;
        }
    }
}


#endif

