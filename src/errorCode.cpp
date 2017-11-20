// Jumping through hoops to remove the dependency
// on the Boost.System library in a way that
// works for all the Boost versions we want to support.
// Using -DBOOST_SYSTEM_NO_DEPRECATED and/or -DBOOST_ERROR_CODE_HEADER_ONLY
// works only with some versions of boost.

#include <boost/system/error_code.hpp>



namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class ErrorCategory : public boost::system::error_category
        {
        public:
          const char * name() const BOOST_SYSTEM_NOEXCEPT {return "ErrorCategoryy";}
          std::string message(int) const BOOST_SYSTEM_NOEXCEPT {return "ErrorCategory";}
        };
        }
}



namespace boost {
    namespace system {
        const error_category & system_category()
        {
          static const ChanZuckerberg::ExpressionMatrix2::ErrorCategory  e;
          return e;
        }

        const error_category & generic_category()
        {
            static const ChanZuckerberg::ExpressionMatrix2::ErrorCategory  e;
            return e;
        }
    }
}
