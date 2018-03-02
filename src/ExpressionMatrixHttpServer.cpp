#include "CellGraph.hpp"
#include "ClusterGraph.hpp"
#include "ExpressionMatrix.hpp"
#include "filesystem.hpp"
#include "randIndex.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
#include "tokenize.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "fstream.hpp"
#include <regex>



// Fill the server function table, which contains the function associated with each keyword.
// The macro can only be used when the keyword and the function name are identical.
#define CZI_ADD_TO_FUNCTION_TABLE(name) serverFunctionTable[string("/") + #name ] = &ExpressionMatrix::name
#define CZI_ADD_TO_FUNCTION_WITH_BROWSER_INFO_TABLE(name) serverFunctionWithBrowserInfoTable[string("/") + #name ] = &ExpressionMatrix::name
void ExpressionMatrix::fillServerFunctionTable()
{
    // Summary.
    serverFunctionTable[""]                                 = &ExpressionMatrix::exploreSummary;
    serverFunctionTable["/"]                                = &ExpressionMatrix::exploreSummary;
    serverFunctionTable["/index"]                           = &ExpressionMatrix::exploreSummary;
    CZI_ADD_TO_FUNCTION_TABLE(exploreHashTableSummary);

    // Genes and gene sets.
    serverFunctionTable["/gene"]                            = &ExpressionMatrix::exploreGene;
    CZI_ADD_TO_FUNCTION_TABLE(compareTwoGenes);
    serverFunctionTable["/geneInformationContent"]          = &ExpressionMatrix::exploreGeneInformationContent;
    serverFunctionTable["/geneSets"]                        = &ExpressionMatrix::exploreGeneSets;
    serverFunctionTable["/geneSet"]                         = &ExpressionMatrix::exploreGeneSet;
    CZI_ADD_TO_FUNCTION_TABLE(removeGeneSet);
    CZI_ADD_TO_FUNCTION_TABLE(createGeneSetFromRegex);
    CZI_ADD_TO_FUNCTION_TABLE(createGeneSetFromGeneNames);
    CZI_ADD_TO_FUNCTION_TABLE(createGeneSetIntersectionOrUnion);
    CZI_ADD_TO_FUNCTION_TABLE(createGeneSetDifference);
    CZI_ADD_TO_FUNCTION_TABLE(createGeneSetUsingInformationContent);

    // Cells and cell sets.
    serverFunctionTable["/cell"]                            = &ExpressionMatrix::exploreCell;
    CZI_ADD_TO_FUNCTION_TABLE(compareTwoCells);
    serverFunctionTable["/cellSets"]                        = &ExpressionMatrix::exploreCellSets;
    serverFunctionTable["/cellSet"]                         = &ExpressionMatrix::exploreCellSet;
    CZI_ADD_TO_FUNCTION_TABLE(createCellSetUsingMetaData);
    CZI_ADD_TO_FUNCTION_TABLE(createCellSetUsingNumericMetaData);
    CZI_ADD_TO_FUNCTION_TABLE(createCellSetIntersectionOrUnion);
    CZI_ADD_TO_FUNCTION_TABLE(createCellSetDifference);
    CZI_ADD_TO_FUNCTION_TABLE(downsampleCellSet);
    CZI_ADD_TO_FUNCTION_TABLE(removeCellSet);

    // Cell meta data.
    serverFunctionTable["/metaData"]                        = &ExpressionMatrix::exploreMetaData;
    CZI_ADD_TO_FUNCTION_TABLE(metaDataHistogram);
    CZI_ADD_TO_FUNCTION_TABLE(metaDataContingencyTable);
    CZI_ADD_TO_FUNCTION_TABLE(removeMetaData);

    // Similar pairs.
    CZI_ADD_TO_FUNCTION_TABLE(similarPairs);
    CZI_ADD_TO_FUNCTION_TABLE(createSimilarPairs);
    CZI_ADD_TO_FUNCTION_TABLE(removeSimilarPairs);

    // Cell graphs.
    serverFunctionTable["/cellGraphs"]                      = &ExpressionMatrix::exploreCellGraphs;
    CZI_ADD_TO_FUNCTION_TABLE(compareCellGraphs);
    serverFunctionTable["/cellGraph"]                       = &ExpressionMatrix::exploreCellGraph;
    CZI_ADD_TO_FUNCTION_TABLE(createCellGraph);
    CZI_ADD_TO_FUNCTION_TABLE(removeCellGraph);

    // Clustering and cluster graphs.
    // CZI_ADD_TO_FUNCTION_TABLE(clusterDialog);
    // CZI_ADD_TO_FUNCTION_TABLE(cluster);
    CZI_ADD_TO_FUNCTION_TABLE(exploreClusterGraphs);
    CZI_ADD_TO_FUNCTION_TABLE(exploreClusterGraph);
    CZI_ADD_TO_FUNCTION_TABLE(exploreClusterGraphSvgWithLabels);
    CZI_ADD_TO_FUNCTION_TABLE(exploreClusterGraphPdfWithLabels);
    nonHtmlKeywords.insert("/exploreClusterGraphPdfWithLabels");
    CZI_ADD_TO_FUNCTION_TABLE(createClusterGraphDialog);
    CZI_ADD_TO_FUNCTION_TABLE(createClusterGraph);
    CZI_ADD_TO_FUNCTION_TABLE(removeClusterGraph);
    CZI_ADD_TO_FUNCTION_TABLE(exploreCluster);
    CZI_ADD_TO_FUNCTION_TABLE(exploreClusterCells);
    CZI_ADD_TO_FUNCTION_TABLE(compareClustersDialog);
    CZI_ADD_TO_FUNCTION_TABLE(compareClusters);
    CZI_ADD_TO_FUNCTION_TABLE(createMetaDataFromClusterGraph);


    // Signature graphs.
    CZI_ADD_TO_FUNCTION_TABLE(exploreSignatureGraphs);
    CZI_ADD_TO_FUNCTION_TABLE(exploreSignatureGraph);
    CZI_ADD_TO_FUNCTION_TABLE(createSignatureGraph);
    CZI_ADD_TO_FUNCTION_TABLE(removeSignatureGraph);


    // Gene graphs.
    CZI_ADD_TO_FUNCTION_TABLE(exploreGeneGraphs);
    CZI_ADD_TO_FUNCTION_WITH_BROWSER_INFO_TABLE(exploreGeneGraph);
    CZI_ADD_TO_FUNCTION_TABLE(createGeneGraph);
    CZI_ADD_TO_FUNCTION_TABLE(removeGeneGraph);
}
#undef CZI_ADD_TO_FUNCTION_TABLE
#undef CZI_ADD_TO_FUNCTION_WITH_BROWSER_INFO_TABLE



// Function that provides simple http functionality
// to facilitate data exploration and debugging.
// It is passed the string of the GET request,
// already parsed using "?=" as separators.
// The request is guaranteed not to be empty.
void ExpressionMatrix::processRequest(
    const vector<string>& request,
    ostream& html,
    const BrowserInformation& browserInformation)
{
    // Look up the keyword to find the function that will process this request.
    const string& keyword = request.front();
    const auto it1 = serverFunctionTable.find(keyword);
    const auto it2 = serverFunctionWithBrowserInfoTable.find(keyword);


    // If we did not find the keyword, see if it is a documentation request.
    if(it1==serverFunctionTable.end() && it2==serverFunctionWithBrowserInfoTable.end()) {



        // See if it is a documentation request of the form help/xyz.
        // Here, xyz is not allowed to contain "../",
        // otherwise we would open a big security hole: requests of the form help/../../file
        // would give access to anywhere in the file system.
        vector<string> tokens;
        boost::algorithm::split(tokens, keyword, boost::algorithm::is_any_of("/"));
        if (keyword.find("../") == string::npos &&
            tokens.front().empty() &&
            tokens[1] == "help") {
            const string name = keyword.substr(6);

            if(serverParameters.docDirectory.empty()) {

                // No documentation directory was specified in the
                // server parameters.
                // satisfy the request from the GitHub Pages.
                // This documentation does not becessarily apply to the
                // version we are running, but it is better than nothing.
                html <<
                    "\r\n"
                    "<!DOCTYPE html>"
                    "<html><script>window.location.replace("
                    "'https://chanzuckerberg.github.io/ExpressionMatrix2/doc/" << name <<
                    "');</script></html>";
                return;


            } else {

                // Satisfy the documentation request as specified
                // in serverParameters.
                ifstream file(serverParameters.docDirectory + "/" + name);
                if (file) {
                    html << "\r\n" << file.rdbuf();
                    return;
                }
            }
        }



        // This is not a documentation request.
        // See if we can interpret as a file name in the current directory.
        // Note that this gives the client access to all files in the
        // current directory.
        if(keyword.size()>1 && keyword[0]=='/') {
            const string fileName = keyword.substr(1);  // Remove the initial slash
            if(fileName.find('/') == string::npos) {    // Make sure there are no other slashes
                ifstream file(fileName);
                if (file) {
                    const string fileExtension = filesystem::extension(fileName);
                    if(fileExtension == "pdf") {
                        html << "Content-Type: application/pdf\r\n";
                    }
                    html << "\r\n" << file.rdbuf();
                    return;
                }
            }
        }



        // We don't know how to satisfy this request.
        // Write a message and stop here.
        html << "\r\nUnsupported request " << keyword;
        return;
    }




    // Write everything that goes before the html body.
    const bool isHtml = nonHtmlKeywords.find(keyword) == nonHtmlKeywords.end();
    if(isHtml) {
        html <<
            "\r\n"
            "<!DOCTYPE html>"
            "<html>"
            "<head>"
            "<link rel=icon href=\"https://s0.wp.com/wp-content/themes/vip/czi/images/build/favicon.ico\" />"
            "<meta charset='UTF-8'>";
        writeStyle(html);
        html <<
            "</head>"
            "<body>";
        writeNavigation(html);
    }



    // We found the keyword. Call the function that processes this keyword.
    // The processing function is only responsible for writing the html body.
    try {
        if(it1 != serverFunctionTable.end()) {
            const auto function = it1->second;
            (this->*function)(request, html);
        } else if(it2 != serverFunctionWithBrowserInfoTable.end()) {
            const auto function = it2->second;
            (this->*function)(request, html, browserInformation);
        }
    } catch(std::exception& e) {
        html << e.what();
    }

    if(isHtml) {
        html << "</body>";
        html << "</html>";
    }
}

ServerParameters::ServerParameters(uint16_t port, string docDirectory) :
    port(port),
    docDirectory(docDirectory)
{
}

void ExpressionMatrix::explore(uint16_t port, const string& docDirectory)
{
    ServerParameters serverParameters(port, docDirectory);
    explore(serverParameters);

}

void ExpressionMatrix::explore(const ServerParameters& serverParameters)
{
    // Store the server parameters.
    this->serverParameters = serverParameters;

    // Validate the docDirectory in the server parameters.
    // If specified and not valid, set it to an empty string.
    if (this->serverParameters.docDirectory.empty()) {
        cout << "A documentation directory was not specified. "
            "Documentation will not be available from the server."  << endl;
    } else {
        ifstream file(serverParameters.docDirectory + "/index.html");
        if (!file) {
            cout << "Documentation index file " << serverParameters.docDirectory;
            cout << "/index.html could not be open. Documentation will not be available from the server." << endl;
            this->serverParameters.docDirectory.clear();
        }
    }

    // Invoke the base class.
    HttpServer::explore(serverParameters.port);
}



void ExpressionMatrix::writeNavigation(ostream& html)
{
    html << "<ul class=navigationMenu>";

    writeNavigation(html, "Genes", {
        {"Genes", "gene"},
        {"Compare two genes", "compareTwoGenes"},
        {"Gene sets", "geneSets"},
        });
    writeNavigation(html, "Cells", {
        {"Cells", "cell"},
        {"Cell meta data", "metaData"},
        {"Compare two cells", "compareTwoCells"},
        {"Cell sets", "cellSets"},
        });
    writeNavigation(html, "Gene graphs", {
        {"Gene graphs", "exploreGeneGraphs"}
        });
    writeNavigation(html, "Cell graphs", {
        {"Find pairs of similar cells", "similarPairs"},
        {"Cell graphs", "cellGraphs"},
        {"Cell clustering", "exploreClusterGraphs"},
        {"Signature graphs", "exploreSignatureGraphs"}
        });
    writeNavigation(html, "Other", {
        {"Run information", "index"},
        {"Hash tables", "exploreHashTableSummary"}
        });

    html << "</ul>";
}


void ExpressionMatrix::writeNavigation(ostream& html, const string& title, const vector<pair <string, string> >& items)
{
    html <<
        "<li class=navigationMenuEntry>"
        "<div class=navigationButton>" << title << "</div>"
        "<div class=navigationItems>";

    for(const auto& item: items) {
        html << "<a class=navigationItem href=" << item.second << ">" << item.first << "</a>";
    }

    html << "</div></li>";

}



#if 0
void ExpressionMatrix::writeNavigation(ostream& html)
{
    writeNavigation(html, "Summary", "index", "");
    writeNavigation(html, "Genes", "gene");
    writeNavigation(html, "Gene sets", "geneSets");
    writeNavigation(html, "Cells", "cell");
    writeNavigation(html, "Compare two cells", "compareTwoCells");
    writeNavigation(html, "Cell sets", "cellSets");
    writeNavigation(html, "Cell meta data", "metaData");
    writeNavigation(html, "Find pairs of similar cells", "similarPairs");
    writeNavigation(html, "Cell graphs", "cellGraphs");
    writeNavigation(html, "Clustering", "exploreClusterGraphs");
    writeNavigation(html, "Signature graphs", "exploreSignatureGraphs");
    writeNavigation(html, "Gene graphs", "exploreGeneGraphs");

    const string helpTooltip =
        serverParameters.docDirectory.empty() ?
        "Click here to read the latest documentation on GitHub "
        "(not necessarily in sync with the version you are using)." :
        "Click here to read the documentation.";
    writeNavigation(html, "Help", "help/index.html", helpTooltip);

    html << "<div style='clear:both;height:10px;'></div>";

}void ExpressionMatrix::writeNavigation(
    ostream& html,
    const string& text,
    const string& url,
    const string& tooltip)
{

    html <<
        "<button onclick=\"window.location = '" << url << "';\"";
    if(!tooltip.empty()) {
        html << " title='" << tooltip << "'";
    }
    html <<
        " onmouseover='this.style.backgroundColor=\"DeepSkyBlue\";'"
        " onmouseout='this.style.backgroundColor=\"LightBlue\";'"
        " style='"
        "background-color:LightBlue;"
        "width:100px;"
        "height:50px;"
        // "box-shadow:0px 0px;"
        // "margin: 0px 1px 0px 0px;"
        // "border-radius:6px;"
        "border-style: none solid solid none;"
        "border-width: 1px;"
        "border-color: White;"
        "text-align:center;"
        "vertical-align:middle;"
        "float:left;"
        "cursor:pointer;"
        "'>\n"
        << text << "</button>";
}
#endif


void ExpressionMatrix::exploreSummary(
    const vector<string>& request,
    ostream& html)
{

    html <<
        "<h1>Expression matrix</h1>"
        "<p>The expression matrix has " << cellCount() <<
        " cells and " << geneCount() <<
        " genes.";


    // If the run directory contains a README.html file, copy it to html.
    if(filesystem::exists("README.html") && filesystem::isRegularFile("README.html")) {
        ifstream readMeFile("README.html");
        // html << "<h1>README.html file</h1>";
        html << readMeFile.rdbuf();
    }



    // If the run directory contains a README file, copy it to html (enclosed in a <pre> element).
    if(filesystem::exists("README") && filesystem::isRegularFile("README")) {
        ifstream readMeFile("README");
        // html << "<h1>README file</h1>
        html << "<pre>";
        html << readMeFile.rdbuf();
        html << "</pre>";
    }



}



void ExpressionMatrix::exploreHashTableSummary(const vector<string>& request, ostream& html)
{
    html << "<h1>Hash table usage analysis</h1>";

    html <<
        "<p>"
        "We use a few hash tables with "
        "<a href='https://en.wikipedia.org/wiki/Open_addressing'>open addressing and linear probing</a>. "
        "These hash tables behave well at low load factors, but their performance "
        "<a href='http://cseweb.ucsd.edu/~kube/cls/100/Lectures/lec16/lec16-27.html'>degrades rapidly</a> "
        "when the load factor is 0.5 or greater. "
        " The table below provides a current usage summary of these hash tables.";

    html <<
        "<table>"
        "<tr><th>Information contained in hash table"
        "<th>Name of controlling creation parameter"
        "<th>Total number of slots"
        "<th>Number of slots in use"
        "<th>Load factor"
        "<th>Performance penalty<br>relative to<br>empty hash table";

    {
        // Hash table used to hold gene names.
        const double capacity = double(geneNames.capacity());
        const double size = double(geneCount());
        const double loadFactor = size / capacity;
        const double penalty = 0.5 * (1. + 1. / ((1.-loadFactor) * (1.-loadFactor))) - 1.;
        const bool isBad = (loadFactor >= 0.3);
        const string color = isBad ? (" style='background-color:red'") : "";
        html <<
            "<tr><td>Gene names"
            "<td class=centered>geneCapacity"
            "<td class=centered>" << capacity <<
            "<td class=centered>" << size;
        const auto oldPrecision = html.precision(3);
        html <<
            "<td class=centered" << color << ">" << loadFactor <<
            "<td class=centered" << color << ">" << penalty;
        html.precision(oldPrecision);
    }

    {
        // Hash table used to hold gene meta data names.
        const double capacity = double(geneMetaDataNames.capacity());
        const double size = double(geneMetaDataNames.size());
        const double loadFactor = size / capacity;
        const double penalty = 0.5 * (1. + 1. / ((1.-loadFactor) * (1.-loadFactor))) - 1.;
        const bool isBad = (loadFactor >= 0.3);
        const string color = isBad ? (" style='background-color:red'") : "";
        html <<
            "<tr><td>Gene meta data names"
            "<td class=centered>geneMetaDataNameCapacity"
            "<td class=centered>" << capacity <<
            "<td class=centered>" << size;
        const auto oldPrecision = html.precision(3);
        html <<
            "<td class=centered" << color << ">" << loadFactor <<
            "<td class=centered>" << penalty;
        html.precision(oldPrecision);
    }

    {
        // Hash table used to hold gene meta data values.
        const double capacity = double(geneMetaDataValues.capacity());
        const double size = double(geneMetaDataValues.size());
        const double loadFactor = size / capacity;
        const double penalty = 0.5 * (1. + 1. / ((1.-loadFactor) * (1.-loadFactor))) - 1.;
        const bool isBad = (loadFactor >= 0.3);
        const string color = isBad ? (" style='background-color:red'") : "";
        html <<
            "<tr><td>Gene meta data values"
            "<td class=centered>geneMetaDataValueCapacity"
            "<td class=centered>" << capacity <<
            "<td class=centered>" << size;
        const auto oldPrecision = html.precision(3);
        html <<
            "<td class=centered" << color << ">" << loadFactor <<
            "<td class=centered>" << penalty;
        html.precision(oldPrecision);
    }

    {
        // Hash table used to hold cell names.
        const double capacity = double(cellNames.capacity());
        const double size = double(cellCount());
        const double loadFactor = size / capacity;
        const double penalty = 0.5 * (1. + 1. / ((1.-loadFactor) * (1.-loadFactor))) - 1.;
        const bool isBad = (loadFactor >= 0.3);
        const string color = isBad ? (" style='background-color:red'") : "";
        html <<
            "<tr><td>Cell names"
            "<td class=centered>cellCapacity"
            "<td class=centered>" << capacity <<
            "<td class=centered>" << size;
        const auto oldPrecision = html.precision(3);
        html <<
            "<td class=centered" << color << ">" << loadFactor <<
            "<td class=centered>" << penalty;
        html.precision(oldPrecision);
    }

    {
        // Hash table used to hold cell meta data names.
        const double capacity = double(cellMetaDataNames.capacity());
        const double size = double(cellMetaDataNames.size());
        const double loadFactor = size / capacity;
        const double penalty = 0.5 * (1. + 1. / ((1.-loadFactor) * (1.-loadFactor))) - 1.;
        const bool isBad = (loadFactor >= 0.3);
        const string color = isBad ? (" style='background-color:red'") : "";
        html <<
            "<tr><td>Cell meta data names"
            "<td class=centered>cellMetaDataNameCapacity"
            "<td class=centered>" << capacity <<
            "<td class=centered>" << size;
        const auto oldPrecision = html.precision(3);
        html <<
            "<td class=centered" << color << ">" << loadFactor <<
            "<td class=centered>" << penalty;
        html.precision(oldPrecision);
    }

    {
        // Hash table used to hold cell meta data values.
        const double capacity = double(cellMetaDataValues.capacity());
        const double size = double(cellMetaDataValues.size());
        const double loadFactor = size / capacity;
        const double penalty = 0.5 * (1. + 1. / ((1.-loadFactor) * (1.-loadFactor))) - 1.;
        const bool isBad = (loadFactor >= 0.3);
        const string color = isBad ? (" style='background-color:red'") : "";
        html <<
            "<tr><td>Cell meta data values"
            "<td class=centered>cellMetaDataValueCapacity"
            "<td class=centered>" << capacity <<
            "<td class=centered>" << size;
        const auto oldPrecision = html.precision(3);
        html <<
            "<td class=centered" << color << ">" << loadFactor <<
            "<td class=centered>" << penalty;
        html.precision(oldPrecision);
    }

    html << "</table>";

}



ostream& ExpressionMatrix::writeMetaDataSelection(
    ostream& html,
    const string& selectName,
    bool multiple) const
{
    set<string> selected;
    return writeMetaDataSelection(html, selectName, selected, multiple);
}
ostream& ExpressionMatrix::writeMetaDataSelection(
    ostream& html,
    const string& selectName,
    const vector<string>& selected,
    bool multiple) const
{
    set<string> selectedSet(selected.begin(), selected.end());
    return writeMetaDataSelection(html, selectName, selectedSet, multiple);
}
ostream& ExpressionMatrix::writeMetaDataSelection(
    ostream& html,
    const string& selectName,
    const set<string>& selected,
    bool multiple) const
{
    // Start the <select> element.
    html << "<select name='" << selectName << "' id='" << selectName << "'";
    if(multiple) {
        html << " multiple";
    }
    html << " style='width:8em;vertical-align:text-top;'>";


    // Write an <option> element for each meta data field.
    CZI_ASSERT(cellMetaDataNamesUsageCount.size() == cellMetaDataNames.strings.size());
    for(StringId i=0; i<cellMetaDataNames.strings.size(); i++) {
        if(cellMetaDataNamesUsageCount[i] == 0) {
            continue;
        }
        const auto cellMetaDataName = cellMetaDataNames.strings[i];
        const string cellMetaDataNameString(cellMetaDataName.begin(), cellMetaDataName.end());
        html << "<option value='" << cellMetaDataNameString << "'";
        if(selected.find(cellMetaDataNameString) != selected.end()) {
            html << " selected=selected";
        }
        html << ">" << cellMetaDataNameString << "</option>";
    }

    // Finish the <select> element.
    html << "</select>";

    return html;
}



void ExpressionMatrix::removeCellGraph(
    const vector<string>& request,
    ostream& html)
{
    string graphName;
    if(!getParameterValue(request, "graphName", graphName)) {
        html << "<p>Missing graph name.";
    } else {
        const auto it = cellGraphs.find(graphName);
        if(it == cellGraphs.end()) {
            html << "<p>Graph " << graphName << " does not exist.";
        } else {
            cellGraphs.erase(it);
            html << "<p>Graph " << graphName << " was removed.";
        }
    }


    html << "<p><form action=cellGraphs><input type=submit autofocus value=Continue></form>";
}



ostream& ExpressionMatrix::writeCellGraphSelection(
    ostream& html,
    const string& selectName,
    bool multiple) const
{
    html << "<select";
    if(multiple) {
        html << " multiple";
    }
    html << " name=" << selectName << "><option value=''></option>";
    for(const auto& p: cellGraphs) {
        const string& graphName = p.first;
        html << "<option value=" << graphName << ">" << graphName << "</option>";
    }
    html << "</select>";

    return html;
}



ostream& ExpressionMatrix::writeNormalizationSelection(
    ostream& html,
    NormalizationMethod selectedMethod) const
    {
    // Begin the select element.
    html << "<select name=normalizationMethod id=normalizationMethod>";

    // Write an <option> element for each possible normalization method.
    for(NormalizationMethod method : validNormalizationMethods) {
        html << "<option value=" << normalizationMethodToShortString(method);
        if(method == selectedMethod) {
            html << " selected=selected";
        }
        html << ">" << normalizationMethodToLongString(method) << "</option>";
    }

    // End the select element.
    html << "</select>";

    return html;
}



NormalizationMethod ExpressionMatrix::getNormalizationMethod(
    const vector<string>& request,
    NormalizationMethod defaultValue)
{
    string normalizationMethodString;
    if(getParameterValue(request, "normalizationMethod", normalizationMethodString)) {
        return normalizationMethodFromShortString(normalizationMethodString);
    } else {
        return defaultValue;
    }
}


#if 0
void ExpressionMatrix::clusterDialog(
    const vector<string>& request,
    ostream& html)
{
    // Get the graph name.
    string graphName;
    if(!getParameterValue(request, "graphName", graphName)) {
        html << "Missing graph name.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }

    // Write the title.
    html << "<h1>Run clustering on graph " << graphName << "</h1>";

    // Form to enter the clustering parameters.
    html <<
        "<form action=cluster>"
        "<h4>Parameters for label propagation clustering</h4>"
        "Random number generator seed: <input type=text name=seed value=231>"
        "<br>Stop after this many iterations without changes: <input type=text name=stableIterationCountThreshold value=3>"
        "<br>Maximum number of iterations: <input type=text name=maxIterationCount value=100>"
        "<br>Meta data name to store the cluster of each cell (leave empty for none): <input type=text name=metaDataName>"
        "<br><h4>Parameters for creation and display of the cluster graph</h4>"
        "In the cluster graph, each vertex corresponds to a cluster of the cell similarity graph."
        "<br>Minimum number of cells for cluster graph vertices: <input type=text name=clusterSizeThreshold value='100'>"
        " Vertices of the cluster graph with fewer than this number of cells will be removed."
        "<br>Similarity threshold for cluster graph edges: <input type=text name=similarityThreshold value='0.5'>"
        " Edges of the cluster graph with similarity lower than this will be removed."
        "<br>Maximum connectivity <input type=text name=maxConnectivity value='3'>"
        " Keep up to this number of best edges for each vertex (this is k of k-NN cluster graph)."
        "<br>Timeout in seconds to compute cluster graph layout <input type=text name=timeout value='30'>"
        "<input type=hidden name=graphName value=" << graphName << ">"
        "<p><input type=submit value='Run clustering' autofocus>"
        "</form>";
}



void ExpressionMatrix::cluster(
    const vector<string>& request,
    ostream& html)
{
    // Get the graph name and the clustering parameters.
    string graphName;
    if(!getParameterValue(request, "graphName", graphName)) {
        html << "Missing graph name.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }
    size_t seed = 231;
    getParameterValue(request, "seed", seed);
    size_t stableIterationCountThreshold = 3;
    getParameterValue(request, "stableIterationCountThreshold", stableIterationCountThreshold);
    size_t maxIterationCount = 50;
    getParameterValue(request, "maxIterationCount", maxIterationCount);
    string metaDataName;
    getParameterValue(request, "metaDataName", metaDataName);
    size_t clusterSizeThreshold;
    getParameterValue(request, "clusterSizeThreshold", clusterSizeThreshold);
    double similarityThreshold;
    getParameterValue(request, "similarityThreshold", similarityThreshold);
    size_t maxConnectivity;
    getParameterValue(request, "maxConnectivity", maxConnectivity);
    double timeout;
    getParameterValue(request, "timeout", timeout);


    // Find the cell graph.
    const auto it = cellGraphs.find(graphName);
    if(it == cellGraphs.end()) {
        html << "<p>Graph " << graphName << " does not exists.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }
    const CellGraphInformation& graphInformation = it->second.first;
    const string& similarPairsName = graphInformation.similarPairsName;
    const SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName);
    const GeneSet& geneSet = similarPairs.getGeneSet();
    CellGraph& graph = *(it->second.second);

    // Write the title.
    html << "<h1>Clustering on graph " << graphName << "</h1>";

    // Do the clustering.
    html << "<pre>";
    graph.labelPropagationClustering(html, seed, stableIterationCountThreshold, maxIterationCount);
    html << "</pre>";


    // If a meta data name was specified, store the  cluster ids in the specified meta data field.
    if(!metaDataName.empty()) {
        storeClusterId(metaDataName, graph);
        html << "<p>The clusters found were stored in cell meta data name " << metaDataName << endl;
    }

    // Create the ClusterGraph.
    ClusterGraph clusterGraph(graph, geneSet);

    // Remove the vertices that correspond to small clusters.
    clusterGraph.removeSmallVertices(clusterSizeThreshold);

    // Compute the average expression for each cluster - that is, for each vertex
    // of the cluster graph.
    BGL_FORALL_VERTICES(v, clusterGraph, ClusterGraph) {
        ClusterGraphVertex& vertex = clusterGraph[v];
        const NormalizationMethod normalizationMethod = NormalizationMethod::L2;    // Use L2 normalization. We might need to make this configurable.
        computeAverageExpression(
            geneSet,
            vertex.cells,
            vertex.averageGeneExpression,
            normalizationMethod);
    }

    // Store in each edge the similarity of the two clusters, computed using the clusters
    // average expression stored in each vertex.
    clusterGraph.computeSimilarities();

    // Remove edges with low similarity.
    clusterGraph.removeWeakEdges(similarityThreshold);    // This may need to be made configurable.
    html << "<p>The cluster graph has " << num_vertices(clusterGraph);
    html << " vertices and " << num_edges(clusterGraph) << " edges.</p>";

    // Make it a k-nn graph.
    clusterGraph.makeKnn(maxConnectivity);

    // Write out the cluster graph in graphviz format.
    clusterGraph.write("ClusterGraph.dot", "", geneNames);



    // Use Graphviz to create a layout of the ClusterGraph in svg format.
    // For better layouts, we would like to use -Goverlap=false because that does not work with the
    // ubuntu graphviz package (it is built without the triangulation library).
    const string command = "timeout " + lexical_cast<string>(timeout) + " sfdp -O -T svg ClusterGraph.dot -Goverlap=scalexy -Gsplines=true";
    const int commandStatus = ::system(command.c_str());
    if(WIFEXITED(commandStatus)) {
        const int exitStatus = WEXITSTATUS(commandStatus);
        if(exitStatus == 124) {
            html <<
                "<p>Computation of the cluster graph layout was interrupted because if was taking too long."
                "<p>You can take one or more of the following steps:<ul>"
                "<li>Decrease the maximum connectivity."
                "<li>Increase the cluster size threshold."
                "<li>Increase the edge similarity threshold."
                "<li>Increase the timeout."
                "</ul>";
            return;
        }
        else if(exitStatus!=0 && exitStatus!=1) {    // sfdp returns 1 all the time just because of the message about missing triangulation.
            throw runtime_error("Error " + lexical_cast<string>(exitStatus) + " running graph layout command: " + command);
        }
    } else if(WIFSIGNALED(commandStatus)) {
        const int signalNumber = WTERMSIG(commandStatus);
        throw runtime_error("Signal " + lexical_cast<string>(signalNumber) + " while running graph layout command: " + command);
    } else {
        throw runtime_error("Abnormal status " + lexical_cast<string>(commandStatus) + " while running graph layout command: " + command);

    }



    // Copy the output svg file to html.
    html <<
        "<h2 id=clusterGraph>Cluster graph</h2>"
        "<p>The cluster graph is shown below. Each vertex represents a cluster. "
        "The most expressed genes in each cluster are listed. "
        "See the table below for details of average gene expression for each cluster.<br>";

    html << ifstream("ClusterGraph.dot.svg").rdbuf();



    // Write to html jQuery and TableSorter so we can make the table below sortable.
    writeJQuery( html);
    writeTableSorter(html);



    // Write a table of average expression for each cluster.
    vector< pair<size_t, ClusterGraph::vertex_descriptor> > sortedVertices;
    BGL_FORALL_VERTICES(v, clusterGraph, ClusterGraph) {
        const ClusterGraphVertex& vertex = clusterGraph[v];
        sortedVertices.push_back(make_pair(vertex.cells.size(), v));
    }
    sort(sortedVertices.begin(), sortedVertices.end(),
        std::greater< pair<size_t, ClusterGraph::vertex_descriptor> >());
    html <<
        "<h2>Average expression in each cluster</h2>"
        "<p>The table contains L2-normalized average expression for each cluster."
        "<p><strong>The table is sortable.</strong> Click on a header to sort by that header. "
        "Click again to reverse the sorting order."
        "<script>"
        "function selectTable() {"
        "    element = document.getElementById(\"averageExpressionTable\");"
        "    selection = window.getSelection();"
        "    range = document.createRange();"
        "    range.selectNodeContents(element);"
        "    selection.removeAllRanges();"
        "    selection.addRange(range);"
        "}"
        "</script>"
        "<p><button onclick='selectTable();'>"
        "Click here to select the entire table.</button>"
        "<p><table id=averageExpressionTable class=tablesorter><thead><tr><th>Gene";
    for(const auto& p: sortedVertices) {
        const ClusterGraph::vertex_descriptor v = p.second;
        const ClusterGraphVertex& vertex = clusterGraph[v];
        html << "<th class=centered>Cluster<br>" << vertex.clusterId << "<br>(" << vertex.cells.size() << " cells)";
    }
    html << "</thead><tbody>";
    for(GeneId localGeneId=0; localGeneId<geneSet.size(); localGeneId++) {
        const GeneId globalGeneId = geneSet.getGlobalGeneId(localGeneId);
        html << "<tr><td class=left>" << geneNames[globalGeneId];
        for(const auto& p: sortedVertices) {
            const ClusterGraph::vertex_descriptor v = p.second;
            const ClusterGraphVertex& vertex = clusterGraph[v];
            html << "<td class=centered>" << vertex.averageGeneExpression[localGeneId];
        }
    }
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#averageExpressionTable').tablesorter();"
        "window.location='#clusterGraph'});"
        "</script>"
        ;


}
#endif



void ExpressionMatrix::createCellGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters.
    string graphName;
    if(!getParameterValue(request, "graphName", graphName)) {
        html << "Missing graph name.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }

    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "Missing cell set name.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }

    string similarPairsName;
    if(!getParameterValue(request, "similarPairsName", similarPairsName)) {
        html << "Missing similar pairs name.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }

    double similarityThreshold;
    if(!getParameterValue(request, "similarityThreshold", similarityThreshold)) {
        html << "Missing or invalid similarity threshold.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }

    int maxConnectivity;
    if(!getParameterValue(request, "maxConnectivity", maxConnectivity)) {
        html << "Missing or invalid max connectivity.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }

    // Check that the name does not already exist.
    if(cellGraphs.find(graphName) != cellGraphs.end()) {
        html << "<p>Graph " << graphName << " already exists.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }


    // Create the graph.
    html << "<div style='font-family:courier'>";
    html << timestamp << "Cell graph creation begins.";
    createCellGraph(graphName, cellSetName, similarPairsName, similarityThreshold, maxConnectivity, false);
    const CellGraphInformation& graphInfo = cellGraphs[graphName].first;
    html <<
        "<br>" << timestamp << "New graph " << graphName << " was created. It has " << graphInfo.vertexCount <<
        " vertices and " << graphInfo.edgeCount << " edges"
        " after " << graphInfo.isolatedRemovedVertexCount << " isolated vertices were removed.";
    html << "</div>";

    // Add a button to continue.
    html << "<p><form action=cellGraphs><input type=submit autofocus value=Continue></form>";


}



// Get a list of the currently available sets of similar pairs.
void ExpressionMatrix::getAvailableSimilarPairs(
    vector<string>& availableSimilarPairs) const
{

    const string fileNamePrefix = directoryName + "/SimilarPairs-";
    const string fileNameSuffix = "-Info";
    const vector<string> directoryContents = filesystem::directoryContents(directoryName);
    for(string name: directoryContents) {
        // Here, name contains the entire file name.
        if(stripPrefixAndSuffix(fileNamePrefix, fileNameSuffix, name)) {
            // Here, now contains just the similar pairs set name.
            availableSimilarPairs.push_back(name);
        }
    }

}


// Get a list of the currently available sets of cell signatures
// ( each corresponding to a possible Lsh objects).
vector<string> ExpressionMatrix::getAvailableLsh() const
{
    vector<string> availableLshNames;

    const string fileNamePrefix = directoryName + "/Lsh-";
    const string fileNameSuffix = "-Info";
    const vector<string> directoryContents = filesystem::directoryContents(directoryName);
    for(string name: directoryContents) {
        // Here, name contains the entire file name.
        if(stripPrefixAndSuffix(fileNamePrefix, fileNameSuffix, name)) {
            // Here, now contains just the similar pairs set name.
            availableLshNames.push_back(name);
        }
    }
    sort(availableLshNames.begin(), availableLshNames.end());

    return availableLshNames;
}


void ExpressionMatrix::exploreMetaData(
    const vector<string>& request,
    ostream& html)
{
    html << "<h2>Histogram a meta data field</h2><form action=metaDataHistogram>";
    html << "Cell set name: ";
    writeCellSetSelection(html, "cellSetName", false);
    html << "<br>Meta data name: ";
    writeMetaDataSelection(html, "metaDataName", set<string>(), false);
    html << "<br><input type=submit value='Create histogram'>";
    html << "</form>";

    html << "<h2>Contingency table for two meta data fields</h2><form action=metaDataContingencyTable>";
    html << "Cell set name: ";
    writeCellSetSelection(html, "cellSetName", false);
    html << "<br>Meta data names: ";
    writeMetaDataSelection(html, "metaDataName0", set<string>(), false);
    html << " and ";
    writeMetaDataSelection(html, "metaDataName1", set<string>(), false);
    html << "<br><input type=submit value='Create contingency table'>";
    html << "</form>";

    html << "<h2>Remove a meta data field</h2><form action=removeMetaData>";
    html << "Cell set name: ";
    writeCellSetSelection(html, "cellSetName", false);
    html << "<br>Meta data name: ";
    writeMetaDataSelection(html, "metaDataName", set<string>(), false);
    html << "<br><input type=submit value='Remove'>";
    html << " Warning: clicking here will remove the specified meta data field from all the cells in the specified cell set. ";
    html << " This operation cannot be undone.";
    html << "</form>";


}



void ExpressionMatrix::metaDataHistogram(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters form the request.
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "<p>Missing cell set name.";
        return;
    }
    string metaDataName;
    if(!getParameterValue(request, "metaDataName", metaDataName)) {
        html << "<p>Missing meta data name name.";
        return;
    }

    // Locate the cell set.
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        html << "<p>Cell set " << cellSetName << " not found.";
        return;
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);

    // Find the StringId corresponding to the specified meta data name.
    const StringId metaDataNameId = cellMetaDataNames(metaDataName);
    if(metaDataNameId == cellMetaDataNames.invalidStringId) {
        html << "<p>Meta data name " << metaDataName << " not found.";
        return;

    }

    // Write a title.
    html << "<h1>Histogram for meta data field " << metaDataName;
    html << " on cell set " << cellSetName << "</h1>";

    // Create the histogram.
    vector< pair<string, size_t> > sortedHistogram;
    histogramMetaData(cellSet, metaDataNameId, sortedHistogram);

    // Write out the histogram.
    html << "<table><tr><th>" << metaDataName << "<th>Frequency";
    for(const auto& p: sortedHistogram) {
        html << "<tr><td>" << p.first << "<td class=centered>" << p.second;
    }
    html << "</table>";
}



void ExpressionMatrix::metaDataContingencyTable(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters form the request.
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "<p>Missing cell set name.";
        return;
    }
    string metaDataName0;
    if(!getParameterValue(request, "metaDataName0", metaDataName0)) {
        html << "<p>Missing meta data name name.";
        return;
    }
    string metaDataName1;
    if(!getParameterValue(request, "metaDataName1", metaDataName1)) {
        html << "<p>Missing meta data name name.";
        return;
    }

    // Locate the cell set.
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        html << "<p>Cell set " << cellSetName << " not found.";
        return;
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);

    // Find the StringId'scorresponding to the specified meta data names.
    const StringId metaDataNameId0 = cellMetaDataNames(metaDataName0);
    if(metaDataNameId0 == cellMetaDataNames.invalidStringId) {
        html << "<p>Meta data name " << metaDataName0 << " not found.";
        return;

    }
    const StringId metaDataNameId1 = cellMetaDataNames(metaDataName1);
    if(metaDataNameId1 == cellMetaDataNames.invalidStringId) {
        html << "<p>Meta data name " << metaDataName1 << " not found.";
        return;

    }

    // Write a title.
    html << "<h1>Contingency table for meta data fields " << metaDataName0;
    html << " and " << metaDataName1;
    html << " on cell set " << cellSetName << "</h1>";

    // Create histograms for the two meta data fields..
    vector< pair<string, size_t> > sortedHistogram0;
    vector< pair<string, size_t> > sortedHistogram1;
    histogramMetaData(cellSet, metaDataNameId0, sortedHistogram0);
    histogramMetaData(cellSet, metaDataNameId1, sortedHistogram1);
    const size_t n0 = sortedHistogram0.size();
    const size_t n1 = sortedHistogram1.size();

    // Create maps that, given a meta data value, give its index in
    // sortedHistogram0 or sortedHistogram1.
    map<string, size_t> map0;
    for(size_t i=0; i<sortedHistogram0.size(); i++) {
        map0.insert(make_pair(sortedHistogram0[i].first, i));
    }
    map<string, size_t> map1;
    for(size_t i=0; i<sortedHistogram1.size(); i++) {
        map1.insert(make_pair(sortedHistogram1[i].first, i));
    }


    // Prepare the contincency table.
    vector< vector<size_t> > matrix(n0, vector<size_t>(n1, 0));

    // Fill in the contingency table.
    for(const CellId cellId: cellSet) {
        const string metaDataValue0 = getCellMetaData(cellId, metaDataNameId0);
        const string metaDataValue1 = getCellMetaData(cellId, metaDataNameId1);
        const size_t i0 = map0[metaDataValue0];
        const size_t i1 = map1[metaDataValue1];
        CZI_ASSERT(i0 < n0);
        CZI_ASSERT(i1 < n1);
        ++(matrix[i0][i1]);
    }


    // Compute the Rand Index and Adjusted Rand Index and write them out.
    double randIndex;
    double adjustedRandIndex;
    computeRandIndex(matrix, randIndex, adjustedRandIndex);
    const auto oldPrecision = html.precision(3);
    html << "<table><tr><th class=left>Rand Index<td class=centered>" << randIndex;
    html << "<tr><th class=left>Adjusted Rand Index<td class=centered>" << adjustedRandIndex << "</table><br>";
    html.precision(oldPrecision);



    // Write out the contingency table.
    html << "<table>";

    // Header row.
    html << "<tr style='vertical-align:top;'><th>" << cellMetaDataNames[metaDataNameId1] << "&#10145;";
    html << "<br>" << cellMetaDataNames[metaDataNameId0] << "&#11015;";
    html << "<th class=left>Total";
    for(size_t i1=0; i1<n1; i1++) {
        html << "<th class=centered>" << sortedHistogram1[i1].first;
    }

    // Row with the totals for metaDataName1.
    html << "<tr><th class=centered>Total";
    html << "<th class=centered>" << cellSet.size();
    for(size_t i1=0; i1<n1; i1++) {
        html << "<th class=centered>" << sortedHistogram1[i1].second;
    }

    // Rows with the contingency table.
    for(size_t i0=0; i0<n0; i0++) {
        const vector<size_t>& row0 = matrix[i0];

        // Name.
        html << "<tr><th class=left>" << sortedHistogram0[i0].first;

        // Total for this metaDanaName0.
        html << "<th class=centered>" << sortedHistogram0[i0].second;

        // Data.
        for(size_t i1=0; i1<n1; i1++) {
            html << "<td class=centered>";
            const size_t frequency = row0[i1];
            if(frequency) {
                html << frequency;
            }
        }
    }
    html << "</table>";
}



void ExpressionMatrix::removeMetaData(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters form the request.
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "<p>Missing cell set name.";
        return;
    }
    string metaDataName;
    if(!getParameterValue(request, "metaDataName", metaDataName)) {
        html << "<p>Missing meta data name name.";
        return;
    }

    // Locate the cell set.
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        html << "<p>Cell set " << cellSetName << " not found.";
        return;
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);

    // Find the StringId corresponding to the specified meta data name.
    const StringId metaDataNameId = cellMetaDataNames(metaDataName);
    if(metaDataNameId == cellMetaDataNames.invalidStringId) {
        html << "<p>Meta data name " << metaDataName << " not found.";
        return;
    }


    // Loop over all cells in the cell set.
    for(const CellId cellId: cellSet) {
        for(auto it=cellMetaData.begin(cellId); it!=cellMetaData.end(cellId); ++it) {
            if((*it).first == metaDataNameId) {
                decrementCellMetaDataNameUsageCount(metaDataNameId);
                cellMetaData.erase(it);
                break;
            }
        }
    }

    html << "<p>Meta data field " << metaDataName;
    html << " was removed from all cells of cell set " << cellSetName;

}
