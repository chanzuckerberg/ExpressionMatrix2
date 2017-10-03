// This file contains member functions of class ExpressionMatrix
// that provide http server functionality related to
// clustering and cluster graphs.

#include "ExpressionMatrix.hpp"
#include "ClusterGraph.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

void ExpressionMatrix::exploreClusterGraphs(
    const vector<string>& request,
    ostream& html)
{
    html <<
        "<h1>Cluster graphs</h1>"
        "<p>A cluster graph describes the result of a clustering operation on a cell graph. "
        "A vertex of the cluster graph corresponds to a cluster. "
        "An edge between two vertices in a cluster graph is created if there is sufficient "
        "similarity between the average expression vectors of the two corresponding clusters. "
        "The similarity value is used to label the edge.";

    if(!clusterGraphs.empty()) {
        html << "<p>The following cluster graphs exist:<ul>";
        for(const auto& p: clusterGraphs) {
            const string& clusterGraphName = p.first;
            html << "<li><a href='exploreClusterGraph?clusterGraphName=" << urlEncode(clusterGraphName) << "'>" << clusterGraphName << "</a>";
        }
        html << "</ul>";
    } else {
        html << "<p>No cluster graphs currently exists.";
        if(cellGraphs.empty()) {
            html <<
                "<p>Before creating a cluster graph you need to create a cell graph. "
                "This can be done <a href=cellGraphs>here</a>.";
            return;
        }
    }

    html <<
        "<p><button onclick=\"window.location = 'createClusterGraphDialog';\">"
        "Run clustering and create a new cluster graph</button>";

}



void ExpressionMatrix::createClusterGraphDialog(
    const vector<string>& request,
    ostream& html)
{
    // Title and explanation.
    html <<
        "<h1>Run clustering and store the result in a cluster graph</h1>"
        "<p>This uses the label propagation algorithm to perform clustering "
        "on an existing cell graph and store the results in a new cluster graph.";

    // Create default-constructed parameters to provide default values in the form below.
    ClusterGraphCreationParameters clusterGraphCreationParameters;


    // Form to enter the clustering parameters.
    html <<
        "<form action=createClusterGraph>"
        "<table>"

        "<tr><th class=left>Cell graph"
        "<td class=centered>";

    writeCellGraphSelection(html, "cellGraphName", false);

    html <<
        "<tr><th class=left>Random number generator seed"
        "<td><input type=text name=seed value='" << clusterGraphCreationParameters.seed << "'>"

        "<tr><th class=left>Stop after this many iterations without changes"
        "<td><input type=text name=stableIterationCount value='" <<
        clusterGraphCreationParameters.stableIterationCount << "'>"

        "<tr><th class=left>Maximum number of iterations"
        "<td><input type=text name=maxIterationCount value='" <<
        clusterGraphCreationParameters.maxIterationCount << "'>"

        "<tr><th class=left>Minimum number of cells for each cluster"
        "<td><input type=text name=minClusterSize value='" <<
        clusterGraphCreationParameters.minClusterSize << "'>"

        "<tr><th class=left>Similarity threshold for cluster graph edges"
        "<td><input type=text name=similarityThreshold value='" <<
        clusterGraphCreationParameters.similarityThreshold << "'>"

        "<tr><th class=left>Maximum connectivity"
        "<td><input type=text name=maxConnectivity value='" <<
        clusterGraphCreationParameters.maxConnectivity << "'>"

        "<tr><th class=left>Timeout in seconds to compute cluster graph layout"
        "<td><input type=text name=timeout value='30'>"

        "<tr><th class=left>Name of the new cluster graph to be created"
        "<td><input type=text name=clusterGraphName required>"

        "</table>"
        "<p><input type=submit value='Run clustering'>"
        "</form>";
}



void ExpressionMatrix::createClusterGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters from the request.
    // Additional parameter validation is done in createClusterGraph.

    string cellGraphName;
    if(!getParameterValue(request, "cellGraphName", cellGraphName)) {
        html << "Missing cell graph name.";
        html << "<p><form action=createClusterGraphDialog><input type=submit value=Continue></form>";
        return;
    }

    ClusterGraphCreationParameters clusterGraphCreationParameters;
    getParameterValue(request, "seed", clusterGraphCreationParameters.seed);
    getParameterValue(request, "stableIterationCount", clusterGraphCreationParameters.stableIterationCount);
    getParameterValue(request, "maxIterationCount", clusterGraphCreationParameters.maxIterationCount);
    getParameterValue(request, "minClusterSize", clusterGraphCreationParameters.minClusterSize);
    getParameterValue(request, "similarityThreshold", clusterGraphCreationParameters.similarityThreshold);
    getParameterValue(request, "maxConnectivity", clusterGraphCreationParameters.maxConnectivity);
    double timeout;
    getParameterValue(request, "timeout", timeout);
    string clusterGraphName;
    if(!getParameterValue(request, "clusterGraphName", clusterGraphName)) {
        html << "Missing cluster graph name.";
        html << "<p><form action=createClusterGraphDialog><input type=submit value=Continue></form>";
        return;
    }

    // Create the cluster graph.
    html << "<pre>";
    createClusterGraph(html, cellGraphName, clusterGraphCreationParameters, clusterGraphName);
    html << "</pre>";

    html << "<p><form action=exploreClusterGraphs><input type=submit value=Continue></form>";
}


void ExpressionMatrix::exploreClusterGraph(
    const vector<string>& request,
    ostream& html)
{
    // Locate the cluster graph.
    string clusterGraphName;
    if(!getParameterValue(request, "clusterGraphName", clusterGraphName)) {
        html << "Missing cluster graph name.";
        html << "<p><form action=exploreClusterGraphs><input type=submit value=Continue></form>";
        return;
    }
    const auto it = clusterGraphs.find(clusterGraphName);
    if(it == clusterGraphs.end()) {
        html << "Cluster graph name " << clusterGraphName << " does not exist.";
        html << "<p><form action=exploreClusterGraphs><input type=submit value=Continue></form>";
        return;
    }
    ClusterGraph& clusterGraph = *(it->second);

    html << "<h1>Cluster graph " << clusterGraphName << "</h1>";
}


