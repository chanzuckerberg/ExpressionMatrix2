#include "Aws.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "boost_lexical_cast.hpp"
#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "iostream.hpp"
#include "stdexcept.hpp"

// Implementation of class AwsS3InputFile.
// See Aws.hpp for more information.

AwsS3InputFile::AwsS3InputFile(const string& path)
{
    // Sanity check on the path.
    if((path.size()<=5) || (path.substr(0, 5) != "s3://")) {
        throw runtime_error("Invalid AWS S3 path " + path);
    }

    // Create the local path using a UUID.
    const string uuid = boost::uuids::to_string(boost::uuids::uuid(boost::uuids::random_generator()()));
    localPath = string("/dev/shm/aws-") + uuid;

    // Copy the file from the AWS path to the local path.
    const string command = "aws s3 cp " + path + " " + localPath;
    const int returnCode = ::system(command.c_str());
    if(returnCode!=0) {
        throw runtime_error("Error " +
            lexical_cast<string>(returnCode) +
            " from command: " + command);
    }

}



AwsS3InputFile::~AwsS3InputFile()
{
    boost::filesystem::remove(localPath);
}
