#include "tokenize.hpp"
#include "fstream.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include "boost_lexical_cast.hpp"

#include "stdexcept.hpp"


// Separate a string into tokens using specified separators.
// This used the boost tokenizer, so it supports quoting, escaping, etc.
// See http://www.boost.org/doc/libs/1_58_0/libs/tokenizer/escaped_list_separator.htm
void ChanZuckerberg::ExpressionMatrix2::tokenize(
    const string& separators,
    const string& inputString,
    vector<string>& tokens,
    bool removeLeadingAndTrailingBlanks)
{
    typedef boost::escaped_list_separator<char> EscapedListSeparator;
    const EscapedListSeparator escapedListSeparator("\\", separators, "\"");

    boost::tokenizer<EscapedListSeparator> tokenizer(inputString, escapedListSeparator);

    tokens.clear();
    tokens.insert(tokens.begin(), tokenizer.begin(), tokenizer.end());

    if(removeLeadingAndTrailingBlanks) {
        for(string& token : tokens) {
            boost::algorithm::trim(token);
        }
    }
}

// Same, but without quoting/escaping.
void ChanZuckerberg::ExpressionMatrix2::tokenizeBare(
    const string& separators,
    const string& inputString,
    vector<string>& tokens)
{
    using Separator = boost::char_separator<char>;
    using Tokenizer = boost::tokenizer<Separator>;

    const Separator separator(separators.c_str());
    Tokenizer tokenizer(inputString, separator);

    tokens.clear();
    tokens.insert(tokens.begin(), tokenizer.begin(), tokenizer.end());
}



// Tokenize all the lines of a file.
void ChanZuckerberg::ExpressionMatrix2::tokenizeFile(
    const string& fileName,
    const string& separators,
    vector<vector<string> >& lines)
{
    // Start with no lines.
    lines.clear();

    // Open the file.
    ifstream file(fileName);

    // Read the entire file.
    string line;
    while(true) {

        // Get a line.
        getline(file, line);

        // If done, stop.
        if(!file) {
            break;
        }

        // Tokenize this line.
        lines.resize(lines.size() + 1);
        tokenize(separators, line, lines.back());
    }
}



// This also checks that  the file is acceptable as an expression count or meta data file:
// - It must have at least two lines.
// - All lines must have the same number of tokens (which cannot be zero), with the possible exception of the
//   first line, which can have one less tokens than all other lines.
void ChanZuckerberg::ExpressionMatrix2::tokenizeFileAndCheck(
    const string& fileName,
    const string& separators,
    vector<vector<string> >& lines)
{
    tokenizeFile(fileName, separators, lines);

    // The file cannot be empty.
    if(lines.empty()) {
        throw runtime_error("File " + fileName + " could not be open or is empty.");
    }

    // The file must have more than one line.
    if(lines.size() == 1) {
        throw runtime_error("File " + fileName + " has only one line. At least two required.");
    }

    // Get the number of tokens of line 1 (line 2, 1-based as reported to the user).
    // It must be at least two.
    const size_t tokenCount = lines[1].size();
    if(tokenCount == 0) {
        throw runtime_error("Line 2 of file " + fileName + " has no tokens. At least one is required.");
    }

    // Check that line 0 (line 1, 1-based as reported to the user)
    // has the same number of tokens or 1 less.
    if(lines[0].size() != tokenCount && lines[0].size() != tokenCount - 1) {
        throw runtime_error("Line 1 of file " + fileName +
            " has " + lexical_cast<string>(lines[0].size()) +
            " tokens, but " + lexical_cast<string>(tokenCount) + "or " +
            lexical_cast<string>(tokenCount - 1) + "  were expected.");
    }

    // Check that all subsequent lines have the same number of tokens as
    // line 1 (line 2, 1-based as reported to the user).
    for(size_t i = 2; i < lines.size(); i++) {
        if(lines[i].size() != tokenCount) {
            throw runtime_error("Line " + lexical_cast<string>(i + 1) + " of file " + fileName +
                " has " + lexical_cast<string>(lines[i].size()) +
                " tokens, but " + lexical_cast<string>(tokenCount) + " were expected.");
        }

    }

}


// Return the number of tokens in the second line of the file with the given name.
// Throw an exception if the file cannot be open or has less than two lines.
size_t ChanZuckerberg::ExpressionMatrix2::countTokensInSecondLine(
    const string& fileName,
    const string& separators)
{
    // Open the file.
    ifstream file(fileName);
    if(!file) {
        throw runtime_error("Could not open file " + fileName);
    }

    // Read the first line.
    string line;
    getline(file, line);
    if(!file) {
        throw runtime_error("Error reading header line of file " + fileName);
    }

    // Read the second line.
    getline(file, line);
    if(!file) {
        throw runtime_error("Error reading the second line of file " + fileName);
    }
    if(line[line.size()-1] == 13) {
        line.resize(line.size()-1); // Remove Windows style line end if necessary.
    }

    // Tokenize it.
    vector<string> tokens;
    tokenize(separators, line, tokens, false);

    // Return the number of tokens.
    return tokens.size();

}



// If the given string is prefix+middle+suffix,
// return true and strip out the prefix and suffix.
// Otherwise, return false and leave the string unchanged.
bool ChanZuckerberg::ExpressionMatrix2::stripPrefixAndSuffix(
    const string& prefix,
    const string& suffix,
    string& s
)
{
    // If it does not begin with the prefix, return false.
    if(!boost::starts_with(s, prefix)) {
        return false;
    }

    // If it does not end with the suffix, return false.
    if(!boost::ends_with(s, suffix)) {
        return false;
    }

    // Strip the prefix and the suffix, then return true.
    s.erase(0, prefix.size());
    s.erase(s.size()-suffix.size(), string::npos);
    return true;

}
