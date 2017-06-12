#ifndef CZI_EXPRESSION_MATRIX2_TOKENIZE_HPP
#define CZI_EXPRESSION_MATRIX2_TOKENIZE_HPP

#include<boost/tokenizer.hpp>

#include "string.hpp"
#include "vector.hpp"


namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        inline void tokenize(
            const string& separators,
            const string& inputString,
            vector<string>& tokens);
    }
}


// See http://www.boost.org/doc/libs/1_58_0/libs/tokenizer/escaped_list_separator.htm
inline void ChanZuckerberg::ExpressionMatrix2::tokenize(
    const string& separators,
    const string& inputString,
    vector<string>& tokens)
{
    typedef boost::escaped_list_separator<char> EscapedListSeparator;
    const EscapedListSeparator escapedListSeparator("\\", separators, "\"");

    boost::tokenizer<EscapedListSeparator> tokenizer(inputString, escapedListSeparator);

    tokens.clear();
    tokens.insert(tokens.begin(), tokenizer.begin(), tokenizer.end());
}


#endif


