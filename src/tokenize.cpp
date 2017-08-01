#include "tokenize.hpp"
#include "fstream.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/algorithm/string/trim.hpp>
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
    	for(string& token: tokens) {
    		boost::algorithm::trim(token);
    	}
    }
}



// Tokenize all the lines of a file.
void ChanZuckerberg::ExpressionMatrix2::tokenizeFile(
	const string& fileName,
    const string& separators,
    vector< vector<string> >& lines)
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
    vector< vector<string> >& lines)
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
	if(lines[0].size()!=tokenCount && lines[0].size()!=tokenCount-1) {
		throw runtime_error("Line 1 of file " + fileName +
			" has " + lexical_cast<string>(lines[0].size()) +
			" tokens, but " + lexical_cast<string>(tokenCount) + "or " +
			lexical_cast<string>(tokenCount-1) + "  were expected.");
	}

	// Check that all subsequent lines have the same number of tokens as
	// line 1 (line 2, 1-based as reported to the user).
	for(size_t i=2; i<lines.size(); i++) {
		if(lines[i].size() != tokenCount) {
			throw runtime_error("Line " + lexical_cast<string>(i+1) + " of file " + fileName +
				" has " + lexical_cast<string>(lines[i].size()) +
				" tokens, but " + lexical_cast<string>(tokenCount) + " were expected.");
		}

	}

}
