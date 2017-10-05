// This file contains member functions of class ExpressionMatrix
// that provide http server functionality related to
// clustering and cluster graphs.

#include "ExpressionMatrix.hpp"
#include "ClusterGraph.hpp"
#include "orderPairs.hpp"
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
            html << "<li><a href='exploreClusterGraph?timeout=30&clusterGraphName=" << urlEncode(clusterGraphName) << "'>" << clusterGraphName << "</a>";
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

    // Get the time out for graph layout computation.
    size_t timeout = 30;
    getParameterValue(request, "timeOut", timeout);

    // Title.
    html << "<h1>Cluster graph " << clusterGraphName << "</h1>";

    // Create the graph layout.
    if(!clusterGraph.computeLayout(timeout, clusterGraphName, geneNames)) {
        throw runtime_error("Error or timeout computing layout for clusger graph " + clusterGraphName);
    }

    // Link to get the output in pdf format.
    html << "<p><a href='exploreClusterGraphPdf?clusterGraphName=" <<
        clusterGraphName <<
        "'>Show this cluster graph in pdf format.</a>";

    // Write the svg to html.
    html << "<p>" << clusterGraph.svg;

}


void ExpressionMatrix::exploreClusterGraphPdf(
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

    // Write out the pdf.
    html << "Content-Type: application/pdf\r\n\r\n" << clusterGraph.pdf;

}



void ExpressionMatrix::exploreCluster(
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

    // Get the cluster id.
    uint32_t clusterId;
    if(!getParameterValue(request, "clusterId", clusterId)) {
        html << "Cluster id is missing.";
        return;
    }

    // Find the vertex of the ClusterGraph corresponding to this cluster.
    const auto jt = clusterGraph.vertexMap.find(clusterId);
    if(jt == clusterGraph.vertexMap.end()) {
        html << "<p>Cluster " << clusterId << " of cluster graph " << clusterGraphName << " does not exist.";
    }
    const ClusterGraph::vertex_descriptor v = jt->second;
    const ClusterGraphVertex& vertex = clusterGraph[v];

    // Title.
    html << "<h1>Cluster " << clusterId << " of cluster graph " << clusterGraphName << "</h1>";
    html << "<p>This cluster has " << vertex.cells.size() << " cells.";

    // Form to display the cells of the cluster with selected meta data.
    html <<
        "<form action=exploreClusterCells>"
        "<input type=submit value='Show the cells of this cluster with the following cell meta data fields:'><br>";
    writeMetaDataSelection(html, "metadata", true);
    html <<
        "<input type=hidden name=clusterGraphName value='" << clusterGraphName << "'>"
        "<input type=hidden name=clusterId value='" << clusterId << "'>"
        "</form>";

    // Write to html jQuery and TableSorter so we can make the table below sortable.
    writeJQuery( html);
    writeTableSorter(html);

    // Sort the average expression for this cluster.
    vector< pair<GeneId, double> > sortedExpression(vertex.averageGeneExpression.size());
    for(GeneId localGeneId=0; localGeneId<vertex.averageGeneExpression.size(); localGeneId++) {
        const GeneId globalGeneId = clusterGraph.geneSet[localGeneId];
        sortedExpression[localGeneId] = make_pair(globalGeneId, vertex.averageGeneExpression[localGeneId]);
    }
    sort(sortedExpression.begin(), sortedExpression.end(), OrderPairsBySecondGreater< pair<GeneId, double> >());

    // Show average expression for this cluster.
    html <<
        "<h2>Average gene expression</h2>"
        "<p>Average expression is computed as the average of the L2-normalized gene "
        "expression vectors for the cells in this cluster. "
        "The average is also L2-normalized."
        "<p>The table only includes genes in the gene set in use for this cluster graph. "
        "All L2 normalizations are also done using only genes in the same gene set."
        "<p><strong>The table below is sortable.</strong> Click on a header to sort by that header. "
        "Click again to reverse the sorting order."
        "<br><table id=expressionTable class=tablesorter><thead>"
        "<tr><th class=centered>Gene<th class=centered>Average<br>expression"
        "</thead><tbody>";
    for(const auto& p: sortedExpression) {
        const GeneId globalGeneId = p.first;
        const double& expression = p.second;
        html << "<tr><td>" << geneNames[globalGeneId] << "<td>" << expression;
    }
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#expressionTable').tablesorter();});"
        "</script>"
        ;

}



void ExpressionMatrix::exploreClusterCells(
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

    // Get the cluster id.
    uint32_t clusterId;
    if(!getParameterValue(request, "clusterId", clusterId)) {
        html << "Cluster id is missing.";
        return;
    }

    // Find the vertex of the ClusterGraph corresponding to this cluster.
    const auto jt = clusterGraph.vertexMap.find(clusterId);
    if(jt == clusterGraph.vertexMap.end()) {
        html << "<p>Cluster " << clusterId << " of cluster graph " << clusterGraphName << " does not exist.";
    }
    const ClusterGraph::vertex_descriptor v = jt->second;
    const ClusterGraphVertex& vertex = clusterGraph[v];



    // Get the names of the meta data to display and the corresponding string ids.
    set<string> metaDataToDisplay;
    getParameterValues(request, string("metadata"), metaDataToDisplay);
    vector< pair<StringId, string> > metaDataToDisplayStrings;
    for(const string& s: metaDataToDisplay) {
        const StringId stringId = cellMetaDataNames(s);
        if(stringId == MemoryMapped::StringTable<StringId>::invalidStringId) {
            html << "<p>Invalid meta data field " << s << " will not be shown.";
        } else {
            metaDataToDisplayStrings.push_back(make_pair(stringId, s));
        }
    }
    // Sort them by string id so they appear in the order in which the meta data
    // was initially created.
    sort(metaDataToDisplayStrings.begin(), metaDataToDisplayStrings.end());



    // Title.
    html << "<h1>Cells of cluster " << clusterId << " of cluster graph " << clusterGraphName << "</h1>";
    html << "<p>This cluster has " << vertex.cells.size() << " cells.";

    // Write to html jQuery and TableSorter so we can make the table below sortable.
    writeJQuery( html);
    writeTableSorter(html);


    // Write out the table with the cells.
    html << "<p><strong>The cell table below is sortable.</strong> Click on a header to sort by that header. "
    "Click again to reverse the sorting order.";
    html << "<br><table id=cellTable class=tablesorter><thead><tr><th class=centered>Cell<br>id<th class=centered>Cell<br>name";
    for(const auto& metaDataFieldName: metaDataToDisplayStrings) {
        html << "<th>" << metaDataFieldName.second;
    }
    html << "</thead><tbody>";
    for(const CellId cellId: vertex.cells) {
        CZI_ASSERT(cellId < cells.size());
        const string& cellName = cellNames[cellId];
        html << "<tr><td class=centered>";
        writeCellLink(html, cellId, true) << "<td class=centered>";
        writeCellLink(html, cellId, false);

        // Write the requested meta data.
        for(const pair<StringId, string>& p: metaDataToDisplayStrings) {
            const StringId metaDateNameStringId = p.first;
            for(const pair<StringId, StringId>& q: cellMetaData[cellId]) {
                if(q.first == metaDateNameStringId) {
                    const StringId metaDataValueStringId = q.second;
                    const auto metaDataValueMemoryRange = cellMetaDataValues(metaDataValueStringId);
                    html << "<td class=centered>";
                    for(const char c: metaDataValueMemoryRange) {
                        html << c;
                    }
                }
            }
        }
    }
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#cellTable').tablesorter();});"
        "</script>"
        ;
}
