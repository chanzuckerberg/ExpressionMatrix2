#ifndef CZI_EXPRESSION_MATRIX2_UUID_HPP
#define CZI_EXPRESSION_MATRIX2_UUID_HPP


#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "string.hpp"
namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        inline string randomUuid()
        {
            using namespace boost::uuids;
            return to_string(uuid(random_generator()()));
        }

    }
}
#endif
