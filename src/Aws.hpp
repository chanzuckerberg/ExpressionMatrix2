// Some code to conveniently deal with files on AWS S3.
#ifndef CZI_EXPRESSION_MATRIX2_AWS_HPP
#define CZI_EXPRESSION_MATRIX2_AWS_HPP

#include "fstream.hpp"
#include "string.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        class AwsS3InputFile;
    }
}


// Class to easily read a file on AWS S3.
// The file is copied locally to /dev/shm, then the local copy is opened for input.
// The local copy is then deleted by the destructor.
// To prevent name collisions, the name of the local copy
// is constructed using a UUID.
// Note that in case of a hard crash the destructor is not called
// and the local copy can be left on /dev/shm.
class ChanZuckerberg::ExpressionMatrix2::AwsS3InputFile : public ifstream {
public:

    // Construct from an AWS S3 path. The path must begin with "s3://",
    // followed by the bucket name.
    // The awscli package must be installed and
    // the user must be configured for access to the bucket via "aws s3" commands.
    // The constructor makes a local copy and opens the local copy.
    AwsS3InputFile(const string& path);

    // The destructor deletes the local copy.
    virtual ~AwsS3InputFile();

private:
    string localPath;
};

#endif
