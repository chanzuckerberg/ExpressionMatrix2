// Definition of macro CZI_ASSERT.
// It is always compiled in, regardless of compilation settings.
// It throws a standard exception if the assertion fails.

#ifndef CZI_EXPRESSION_MATRIX2_CZI_ASSERT_H
#define CZI_EXPRESSION_MATRIX2_CZI_ASSERT_H

#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <string>

// Gcc (for backtraces).
#include "execinfo.h"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        inline void writeBackTrace();
    }
}

#define CZI_ASSERT(expression) ((expression) ? (static_cast<void>(0)) : \
    (/*writeBackTrace(),*/ throw std::runtime_error(std::string("Assertion failed: ") + #expression + " at " + BOOST_CURRENT_FUNCTION + " in " +  __FILE__ + " line " + boost::lexical_cast<std::string>(__LINE__))))


#if 0
inline void ChanZuckerberg::ExpressionMatrix2::writeBackTrace()
{
    const int bufferSize = 64;  // To avoid extremely long, useless backtraces.
    void* buffer[bufferSize];
    ::backtrace(buffer, bufferSize);
    ::backtrace_symbols_fd(buffer, bufferSize, ::fileno(::stdout));
}
#endif


#endif

