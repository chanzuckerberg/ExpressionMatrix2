// To eliminate the runtime dependency on the boost system library,
// I pasted here the contents of boost_1_58_0/libs/system/error_code.cpp
// Alternative solutions suggested in the web
// consist of defining BOOST_SYSTEM_NO_DEPRECATED or
// BOOST_ERROR_CODE_HEADER_ONLY. Neither of those worked.



//  error_code support implementation file  ----------------------------------//

//  Copyright Beman Dawes 2002, 2006

//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See library home page at http://www.boost.org/libs/system

//----------------------------------------------------------------------------//

// define BOOST_SYSTEM_SOURCE so that <boost/system/config.hpp> knows
// the library is being built (possibly exporting rather than importing code)
#define BOOST_SYSTEM_SOURCE

#include <boost/system/error_code.hpp>

#ifndef BOOST_ERROR_CODE_HEADER_ONLY
#include <boost/system/detail/error_code.ipp>
#endif
