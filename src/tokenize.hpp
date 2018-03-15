#ifndef CZI_EXPRESSION_MATRIX2_TOKENIZE_HPP
#define CZI_EXPRESSION_MATRIX2_TOKENIZE_HPP

#include<boost/tokenizer.hpp>

#include "string.hpp"
#include "vector.hpp"


namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        // With quoting and escaping.
        void tokenize(
            const string& separators,
            const string& inputString,
            vector<string>& tokens,
            bool removeLeadingAndTrailingBlanks = false);

        // Without quoting and escaping.
        void tokenizeBare(
            const string& separators,
            const string& inputString,
            vector<string>& tokens);

        void tokenizeFile(
            const string& fileName,
            const string& separators,
            vector<vector<string> >& lines);



        // This also checks that  the file is acceptable as an expression count or meta data file:
        // - It must have at least two lines.
        // - All lines must have the same number of tokens (which cannot be zero), with the possible exception of the
        //   first line, which can have one less tokens than all other lines.
        void tokenizeFileAndCheck(
            const string& fileName,
            const string& separators,
            vector<vector<string> >& lines);

        // Return the number of tokens in the second line of the file with the given name.
        // Throw an exception if the file cannot be open or has less than two lines.
        size_t countTokensInSecondLine(
            const string& fileName,
            const string& separators);

        inline void removeWindowsLineEnd(string& s)
        {
            if(!s.empty() && s.back()==13) {
                s.resize(s.size()-1);
            }
        }



        // If the given string is prefix+middle+suffix,
        // return true and strip out the prefix and suffix.
        // Otherwise, return false and leave the string unchanged.
        bool stripPrefixAndSuffix(
            const string& prefix,
            const string& suffix,
            string&
        );

    }
}

#endif


