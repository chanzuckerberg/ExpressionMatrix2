#include "ExpressionMatrix.hpp"
#include "string.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



void writeHelp()
{
    cout <<
        "Usage:\n\n"
        "--help (or -h)\n"
        "Write this help.\n\n"
        "--port portNumber (or -p portNumber)\n"
        "First port number to try (default 17100).\n\n"
        "--data dataDirectory (or -d dataDirectory)\n"
        "Directory containing run data (default \"data\").\n"
        "Will be created if it does not exist.\n\n"
        "--allow-remote (or -r)\n"
        "Specify this to allow remote connections.\n"
        "If not specified, only local connections are accepted\n"
        "(the browser must run on the same machine).\n\n"
        "Invoking without arguments is equivalent to \"--port 17100 --data data\".\n\n"
        "While this is running, start a browser on the same machine\n"
        "and point it to this URL: \"localhost:17100\",\n"
        "replacing the port number if necessary.\n";
}

int main(int argumentCount, const char** arguments)
{

    // Process the arguments.
    uint16_t portNumber = 17100;
    string directory = "data";
    bool localOnly = true;
    int i = 1;
    while(i<argumentCount) {
        using std::string;
        const char* argument = arguments[i++];



        // Process "--help".
        if(
            argument == string("--help") ||
            argument == string("-h")
            ) {
            writeHelp();
            return 1;
        }



        // Process "--port".
        if(
            argument == string("--port") ||
            argument == string("-p")
            ) {
            if(i == argumentCount) {
                cout << "Port number is missing." << endl;
                return 1;
            }
            try {
                portNumber = boost::lexical_cast<uint16_t>(arguments[i++]);
            } catch(boost::bad_lexical_cast) {
                cout << "Invalid port number." << endl;
                return 1;
            }
            continue;
        }



        // Process "--data".
        if(
            argument == string("--data") ||
            argument == string("-d")
            ) {
            if(i == argumentCount) {
                cout << "Data directory is missing." << endl;
                return 1;
            }
            directory = arguments[i++];
            continue;
        }



        // Process "--allow-remote".
        if(
            argument == string("--allow-remote") ||
            argument == string("-r")
            ) {
            localOnly = false;
            continue;
        }



        // If we got here, we have an option we don't understand.
        cout << "Invalid argument or option: " << argument << endl;
        writeHelp();
        return 1;
    }


    // Run.
    try {
        ExpressionMatrix e(directory);
        e.explore(portNumber, "", localOnly);
    } catch(runtime_error& e) {
        cout << e.what();
        return 0;
    }
    
    return 0;
}
