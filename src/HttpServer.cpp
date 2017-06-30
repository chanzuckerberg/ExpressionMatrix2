// Implementation of class HttpServer - see HttpServer.hpp for more information.

#include "HttpServer.hpp"
#include "sstream.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/algorithm/string.hpp>

#include <boost/asio/io_service.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ip/v6_only.hpp>
using namespace boost;
using namespace asio;
using namespace ip;

#include "iostream.hpp"
#include <sstream>



// This function puts the server into an endless loop
// of processing requests.
// This is trhe function that the base class should call to start the server.
void HttpServer::explore(uint16_t port)
{
    // Create the acceptor, making sure to accept both ipv4 and ipv6 ip addresses.
    io_service service;
    tcp::acceptor acceptor(service);
    tcp::endpoint endpoint(tcp::v6(), port);
    acceptor.open(endpoint.protocol());
    v6_only ipv6Option(false);
    acceptor.set_option(ipv6Option);
    acceptor.bind(endpoint);
    acceptor.listen();
    cout << "Listening for http requests on port " << port << endl;

    // Endless loop over incoming connections.
    while(true) {
          tcp::iostream s;
          tcp::endpoint remoteEndpoint;
          boost::system::error_code errorCode;
          acceptor.accept(*s.rdbuf(), remoteEndpoint, errorCode);
          if(errorCode) {
              // If interrupted with Ctrl-C, we get here.
              cout << "\nError code from accept: " << errorCode.message() << endl;
              s.close();        // Should not be necessary.
              acceptor.close(); // Should not be necessary
              return;
          }
          cout << timestamp << remoteEndpoint.address().to_string() << " ";
          processRequest(s);
    }
}



void HttpServer::processRequest(tcp::iostream& s)
{
    // Get the first line, which must contain the GET request.
    string requestLine;
    getline(s, requestLine);

    // Parse it to get only the request string portion.
    // It is the second word of the first line.
    vector<string> tokens;
    boost::algorithm::split(tokens, requestLine, boost::algorithm::is_any_of(" "));
    if(tokens.size() != 3) {
        s << "Unexpected number of tokens in http request: expected 3, got " << tokens.size();
        cout << "Unexpected number of tokens in http request: expected 3, got " << tokens.size();
        return;
    }
    if(tokens.front() != "GET") {
        s << "Unexpected keyword in http request: " << tokens.front();
        cout << "Unexpected keyword in http request: " << tokens.front();
        return;
    }
    const string& request = tokens[1];
    if(request.empty()) {
        s << "Empty GET request: " << requestLine;
        cout << "Empty GET request: " << requestLine;
        return;
    }

    // Parse the request.
    cout << requestLine << endl;
    boost::algorithm::split(tokens, request, boost::algorithm::is_any_of("?=&"));

    // Do URl decoding on each token.
    // This takes care of % encoding, which the browser will do if it has to send special characters.
    // Note that we have to do this after parsing the request into tokens.
    // With this, we can support special characters in cell meta data, cell set names, graph names, etc.
    for(string& token: tokens) {
    	string newToken;
    	urlDecode(token, newToken);
    	if(newToken != token) {
    		cout << "Request token " << token << " decoded as " << newToken << endl;
    	}
    	token = newToken;
    }

    // Read the rest of the input from the client, but ignore it.
    // We have to do this, otherwise the client may get a timeout.
    string line;
    while(true) {
        if(!s) {
            break;
        }
        getline(s, line);
        if(!s) {
            break;
        }
        // cout << "Got another line: " << line.size() /* << "<<<" << line << ">>>" */ << endl;
        if(line.size()==1) {
            break;
        }
    }

    // Write the success response.
    // We don't write the required empty line, so the derived class can send headers
    // if it wants to.
    s << "HTTP/1.1 200 OK\r\n";

    // The derived class processes the request.
    processRequest(tokens, s);
}



// Return all values assigned to a parameter.
// For example, if the request has ...&a=xyz&a=uv,
// when called with argument "a" returns a set containing "xyz" and "uv".
void HttpServer::getParameterValues(const vector<string>& request, const string& name, set<string>& values)
{
    for(size_t i=0; i<request.size()-1; i++) {
        if(request[i]==name) {
            values.insert(request[i+1]);
        }
    }

}
void HttpServer::getParameterValues(const vector<string>& request, const string& name, vector<string>& values)
{
    for(size_t i=0; i<request.size()-1; i++) {
        if(request[i]==name) {
            values.push_back(request[i+1]);
        }
    }

}


// This takes care of percent encoding in the url.
// I adapted it from boost asio example http/server3.
bool HttpServer::urlDecode(const string& in, string& out)
{
  out.clear();
  out.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i)
  {
    if (in[i] == '%')
    {
      if (i + 3 <= in.size())
      {
        int value = 0;
        std::istringstream is(in.substr(i + 1, 2));
        if (is >> hex >> value)
        {
          out += static_cast<char>(value);
          i += 2;
        }
        else
        {
          return false;
        }
      }
      else
      {
        return false;
      }
    }
    else if (in[i] == '+')
    {
      out += ' ';
    }
    else
    {
      out += in[i];
    }
  }
  return true;
}



// Function to do URL encoding.
// This is necessary when we create a hyperlink that contains special characters.
// Adapted from here:
// https://stackoverflow.com/questions/154536/encode-decode-urls-in-c
string HttpServer::urlEncode(const string& s)
{
    ostringstream escaped;
    escaped.fill('0');
    escaped << hex;

    for (string::const_iterator i = s.begin(), n = s.end(); i != n; ++i) {
        string::value_type c = (*i);

        // Keep alphanumeric and other accepted characters intact
        if (isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~') {
            escaped << c;
            continue;
        }

        // Any other characters are percent-encoded
        escaped << std::uppercase;
        escaped << '%' << std::setw(2) << int((unsigned char) c);
        escaped << std::nouppercase;
    }

    return escaped.str();
}



void HttpServer::writeStyle(ostream& html)
{
    html << R"%(
<style>
    * {
        font-family: sans-serif;
    }
    pre {
        font-family: courier;
    }
    p, input {
        font-size: 16px;
    }
    h1, h2, h3 {
        color: LightCoral;
    }
    table {
        border-collapse: collapse;
    }
    th, td {
        border: 2px solid LightCoral;
    }
    th {
        font-weight: bold;
        text-align: center;
    }
    th.left {
        text-align: left;
    }
    td.centered {
        text-align: center;
    }
    a {
        color: LightCoral;
    }
</style>
    )%";
}



ostream& HttpServer::writeJQuery(ostream& html)
{
    html << "<script src='http://ajax.googleapis.com/ajax/libs/jquery/1.8.3/jquery.min.js'></script>";
    return html;

}



ostream& HttpServer::writeTableSorter(ostream& html)
{
    html << "<script src='http://tablesorter.com/__jquery.tablesorter.min.js'></script>";
    return html;
}


