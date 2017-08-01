#ifndef CZI_EXPRESSION_MATRIX2_TOKENIZE_HPP
#define CZI_EXPRESSION_MATRIX2_TOKENIZE_HPP

#include<boost/tokenizer.hpp>

#include "string.hpp"
#include "vector.hpp"


namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

    void tokenize(
        const string& separators,
        const string& inputString,
        vector<string>& tokens,
		bool removeLeadingAndTrailingBlanks = false);

    void tokenizeFile(
    	const string& fileName,
        const string& separators,
        vector< vector<string> >& lines);



    // This also checks that  the file is acceptable as an expression count or meta data file:
    // - It must have at least two lines.
    // - All lines must have the same number of tokens (which cannot be zero), with the possible exception of the
    //   first line, which can have one less tokens than all other lines.
    void tokenizeFileAndCheck(
    	const string& fileName,
        const string& separators,
        vector< vector<string> >& lines);

    }
}

#endif


