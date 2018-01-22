// This file contains member functions of class ExpressionMatrix
// that provide http server functionality related to
// clustering and cluster graphs.

#include "ExpressionMatrix.hpp"
#include "ClusterGraph.hpp"
#include "orderPairs.hpp"
#include "tokenize.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "iterator.hpp"

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
        html << "<p>The following cluster graphs exist:<table>";
        for(const auto& p: clusterGraphs) {
            const string& clusterGraphName = p.first;
            html <<
                "<tr><td><a href='exploreClusterGraph?timeout=30&clusterGraphName=" <<
                urlEncode(clusterGraphName) << "'>" << clusterGraphName << "</a>"
                "<td class=centered><form action=removeClusterGraph>"
                "<input type=text hidden name=clusterGraphName value='" << clusterGraphName <<
                "'><input type=submit value=Remove></form>";
        }
        html << "</table>";
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

        "<tr><th class=left>Similarity threshold to remove cluster graph edges"
        "<td><input type=text name=similarityThreshold value='" <<
        clusterGraphCreationParameters.similarityThreshold << "'>"

        "<tr><th class=left>Similarity threshold to merge cluster graph vertices"
        "<td><input type=text name=similarityThresholdForMerge value='" <<
        clusterGraphCreationParameters.similarityThresholdForMerge << "'>"

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
    getParameterValue(request, "similarityThresholdForMerge", clusterGraphCreationParameters.similarityThresholdForMerge);
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

    // Links to get the layout with labels in svg or pdf format.
    html << "<p>Show this cluster graph with vertex labels in "
        "<a href='exploreClusterGraphSvgWithLabels?clusterGraphName=" <<
        clusterGraphName <<
        "'>svg</a>";
    html << " or <a href='exploreClusterGraphPdfWithLabels?clusterGraphName=" <<
        clusterGraphName <<
        "'>pdf</a> format.";

    // Link to compare gene expression in sets of clusters
    html << "<br><a href='compareClustersDialog?clusterGraphName=" <<
        clusterGraphName <<
        "'>Compare gene expression between clusters.</a>";

    // Link to create meta data from the cluster ids in this cluster graph.
    html <<
        "<form action=createMetaDataFromClusterGraph>"
        "<input type=submit value='Store the cluster ids in this cluster graph in cell meta data name'> "
        "<input type=text name=metaDataName>"
        "<input type=hidden name=clusterGraphName value=" << clusterGraphName << ">"
        "</form>";

    // Write the svg layout without labels to html,
    // computing it first if necessary.
    clusterGraph.computeLayout(timeout, clusterGraphName, geneNames, false);
    html << "<p>" << clusterGraph.svgLayoutWithoutLabels;

}





void ExpressionMatrix::removeClusterGraph(
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
        html << "Cluster graph " << clusterGraphName << " does not exist.";
        html << "<p><form action=exploreClusterGraphs><input type=submit value=Continue></form>";
        return;
    }

    clusterGraphs.erase(it);
    html << "Cluster graph " << clusterGraphName << " was removed.";
    html << "<p><form action=exploreClusterGraphs><input type=submit value=Continue></form>";

}



void ExpressionMatrix::exploreClusterGraphPdfWithLabels(
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

    // Compute the layouts with labels, if needed.
    const int timeoutSeconds = 30;
    clusterGraph.computeLayout(timeoutSeconds, clusterGraphName, geneNames, true);

    // Write out the pdf layout with labels.
    html << "Content-Type: application/pdf\r\n\r\n" << clusterGraph.pdfLayoutWithLabels;

}



void ExpressionMatrix::exploreClusterGraphSvgWithLabels(
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

    // Compute the layouts with labels, if needed.
    const int timeoutSeconds = 30;
    clusterGraph.computeLayout(timeoutSeconds, clusterGraphName, geneNames, true);

    // Write out the svg layout with labels.
    html << "<h1>Cluster graph " << clusterGraphName << "</h1>";
    html << clusterGraph.svgLayoutWithLabels;

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
        "<h2>Show the cells of this cluster</h2>"
        "<form action=exploreClusterCells>"
        "<input type=submit value='Show the cells of this cluster'>"
        "<br>Select zero, one, or multiple cell meta data fields to be shown:<br>";
    writeMetaDataSelection(html, "metadata", true);
    html <<
        "<br>Also include expression for the following genes:"
        "<br><input type=text name=genes size=40>"
        "<br>(enter gene names separated by spaces)."
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

    // Get the genes for which we will write expression counts.
    string genesString;
    getParameterValue(request, "genes", genesString);
    vector<string> requestedGeneNames;
    tokenize(" ", genesString, requestedGeneNames, true);

    // Find the corresponding gene ids.
    vector<GeneId> geneIds;
    for(const string& geneName: requestedGeneNames) {
        const GeneId geneId = geneIdFromName(geneName);
        if(geneId == invalidGeneId) {
            html << "<p>Expression counts for non-existent gene " << geneName << " will not be shown.";
        } else {
            geneIds.push_back(geneId);
        }
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
    for(const GeneId geneId: geneIds) {
        html << "<th>" << geneNames[geneId];
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

        // Write the requested expression counts.
        for(const GeneId geneId: geneIds) {
            html << "<td class=centered>" << getCellExpressionCount(cellId, geneId);
        }
    }
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#cellTable').tablesorter();});"
        "</script>"
        ;
}



void ExpressionMatrix::compareClustersDialog(
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


    html <<
        "<h1>Compare cluster gene expression in cluster graph " << clusterGraphName << "</h1>"
        "<p>Multiple cluster ids can be entered separated by space."
        "<form action=compareClusters><table>"

        "<tr><th class=left>Clusters in the first set to compare"
        "<td>Specify clusters: <input type=text name=aClusters>"

        "<tr><th class=left rowspan=3>Clusters in the second set to compare"
        "<td class=left><input type=radio name=bClustersOption value=specifiedClusters checked=checked>Specify clusters:"
        "<input type=text name=bClusters>"
        "<tr><td class=left><input type=radio name=bClustersOption value=allOther>All clusters not in the first set"
        "<tr><td class=left><input type=radio name=bClustersOption value=entireGraph>Entire cluster graph"

        "</table>"
        "<input type=checkbox name=detailed>Write detailed tabular output"
        "<input type=hidden name=clusterGraphName value='" << clusterGraphName << "'>"
        "<br><input type=submit value=Compare>"
        "</form>";
}



void ExpressionMatrix::createMetaDataFromClusterGraph(
    const vector<string>& request,
    ostream& html)
{
    // Locate the cluster graph.
    string clusterGraphName;
    getParameterValue(request, "clusterGraphName", clusterGraphName);

    // Get the name of the meta data field where the cluster ids should be stored.
    string metaDataName;
    getParameterValue(request, "metaDataName", metaDataName);
    if(metaDataName.empty()) {
        throw runtime_error("Missing meta data name.");
    }

    createMetaDataFromClusterGraph(clusterGraphName, metaDataName);
    html << "<p>The cluster ids in cluster graph " << clusterGraphName <<
        " were stored in cell meta data " << metaDataName << ".";
}



void ExpressionMatrix::compareClusters(
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



    // Get the clusters in the first set.
    string aClustersString;
    getParameterValue(request, "aClusters", aClustersString);
    vector<string> tokens;
    tokenize(" ", aClustersString, tokens, true);
    array< vector<uint32_t>, 2> clusterIds;
    vector<uint32_t>& aClusterIds = clusterIds[0];
    for(const string& token: tokens) {
        try {
            const uint32_t clusterId = lexical_cast<uint32_t>(token);
            if(clusterGraph.vertexMap.find(clusterId) == clusterGraph.vertexMap.end()) {
                html << "<p>Invalid cluster id " << token << " ignored.";
            } else {
                aClusterIds.push_back(lexical_cast<uint32_t>(token));
            }
        } catch(bad_lexical_cast) {
            html << "<p>Invalid cluster id " << token << " ignored.";
        }
    }



    // Get the clusters in the second set.
    string bClustersOption = "specifiedClusters";
    getParameterValue(request, "bClustersOption", bClustersOption);
    vector<uint32_t>& bClusterIds = clusterIds[1];
    if(bClustersOption == "specifiedClusters") {
        string bClustersString;
        getParameterValue(request, "bClusters", bClustersString);
        vector<string> tokens;
        tokenize(" ", bClustersString, tokens, true);
        for(const string& token: tokens) {
            try {
                const uint32_t clusterId = lexical_cast<uint32_t>(token);
                if(clusterGraph.vertexMap.find(clusterId) == clusterGraph.vertexMap.end()) {
                    html << "<p>Invalid cluster id " << token << " ignored.";
                } else {
                    bClusterIds.push_back(lexical_cast<uint32_t>(token));
                }
            } catch(bad_lexical_cast) {
                html << "<p>Invalid cluster id " << token << " ignored.";
            }
        }
    }
    if(bClustersOption == "allOther") {
        for(const auto& p: clusterGraph.vertexMap) {
            const uint32_t clusterId = p.first;
            if(find(aClusterIds.begin(), aClusterIds.end(), clusterId) ==aClusterIds.end()) {
                bClusterIds.push_back(p.first);
            }
        }
    }
    if(bClustersOption == "entireGraph") {
        for(const auto& p: clusterGraph.vertexMap) {
            bClusterIds.push_back(p.first);
        }
    }

    // See if detailed tabular output was requested.
    string detailOption = "specifiedClusters";
    getParameterValue(request, "detailed", detailOption);


    // Write a title.
    html << "<h1>Compare cluster gene expression in cluster graph " << clusterGraphName << "</h1>";

    // Write the clusters in the sets being compared.
    html << "<p>Clusters in first set being compared: ";
    copy(aClusterIds.begin(), aClusterIds.end(), ostream_iterator<uint32_t>(html, " "));
    html << "<br>Clusters in second set being compared: ";
    if(bClustersOption == "allOther") {
        html << "all clusters not in the first set, ";
    }
    if(bClustersOption == "entireGraph") {
        html << "entire cluster graph, ";
    }
    copy(bClusterIds.begin(), bClusterIds.end(), ostream_iterator<uint32_t>(html, " "));

    if(aClusterIds.empty() || bClusterIds.empty()) {
        html << "<p>The two sets of cluster cannot be empty.";
        return;
    }



    // Compute average expression for the two sets of clusters.
    // The averages are weighted by the number of cells in each cluster,
    // then L2-normalized.
    array< vector<double>, 2> setsAverageExpression;
    const size_t geneCount = clusterGraph.geneSet.size();
    for(int i=0; i<2; i++) {

        // Initialize the average expression for this set.
        setsAverageExpression[i].resize(geneCount, 0.);

        // Loop over the clusters in this set.
        for(uint32_t clusterId: clusterIds[i]) {
            const ClusterGraph::vertex_descriptor v = clusterGraph.vertexMap[clusterId];
            const ClusterGraphVertex& vertex = clusterGraph[v];
            const double weight = double(vertex.cells.size());
            for(size_t j=0; j<geneCount; j++) {
                setsAverageExpression[i][j] += weight * vertex.averageGeneExpression[j];
            }
        }

        // L2-normalize the average gene expression for this set.
        double sum2 = 0.;
        for(const double x: setsAverageExpression[i]) {
            sum2+= x*x;
        }
        const double factor = 1./sqrt(sum2);
        for(double& x: setsAverageExpression[i]) {
            x *= factor;
        }
    }



    // Write a scatter plot of the expressions of the two sets of clusters.
    html <<
        "<br>In the scatter plot below, expression values less than 10<sup>-6</sup> are replaced with 10<sup>-6</sup>."
        "<br>All expression values are L2-normalized."
        "<br>You can hover on a data point to display the gene name."
        "<script src='https://www.gstatic.com/charts/loader.js'></script>"
        "<script>"
        "    google.charts.load('current', {'packages':['corechart']});"
        "    google.charts.setOnLoadCallback(drawChart);"

        "    var data;"
        "    var chart;"
        "    var options = {"
        "        chartArea: {top:10, right:10, bottom:100, left:100},"
        "        hAxis: {title: 'First set', logScale: true, minValue:0.000001, maxValue:1},"
        "        vAxis: {title: 'Second set', logScale: true, minValue:0.000001, maxValue:1},"
        "        legend: 'none',"
        "        pointSize: 1,"
        "        tooltip: {isHtml: true, trigger: 'both'}"
        "    };"


        "    function drawChart()"
        "    {"
        "        data = new google.visualization.DataTable();"
        "        data.addColumn('number', 'x');"
        "        data.addColumn('number', 'y');"
        "        data.addColumn({type: 'string', role: 'tooltip', 'p': {'html': true}});"
        "        data.addRows([";
    for(size_t j=0; j<geneCount; j++) {
        const GeneId globalGeneId = clusterGraph.geneSet[j];
        const double count0 = max(1.e-6, setsAverageExpression[0][j]);
        const double count1 = max(1.e-6, setsAverageExpression[1][j]);
        html << "[" << count0 << "," << count1 << ",";
        html << "\"<a href='gene?geneId=" << globalGeneId << "'>" << geneNames[globalGeneId] << "</a>\"";
        html << "],";
    }
    html <<
        "        ]);"
        "        chart = new google.visualization.ScatterChart(document.getElementById('scatterPlot'));"
        "        chart.draw(data, options);"
        "    }"

        "</script>"
        "<div id='scatterPlot' style='float: left;width: 600px; height: 600px;'></div>"
        "<div style='float:left'>"
        "</div>"
        "<div style='clear:both;' />"
        ;



    // If detailed output was not requested, stop here.
    if(detailOption != "on") {
        return;
    }

    // Write to html jQuery and TableSorter so we can make the table below sortable.
    writeJQuery( html);
    writeTableSorter(html);

    // Show the average expression of the two sets of clusters in a table.
    html <<
        "<h2>Average gene expression</h2>"
        "<p><strong>The table below is sortable.</strong> Click on a header to sort by that header. "
        "Click again to reverse the sorting order."
        "<br><table id=expressionTable class=tablesorter><thead><tr>"
        "<th class=centered>Gene"
        "<th class=centered>First set"
        "<th class=centered>Second set"
        "<th class=centered>Max<br>(first set,<br>second set)"
        "<th class=centered>Ratio<br>first set /<br>second set"
        "<th class=centered>Ratio<br>second set /<br>first set"
        "</thead><tbody>";
    // const auto oldPrecision = html.precision(6);
    const auto oldOptions = html.setf(std::ios::fixed);
    for(size_t j=0; j<geneCount; j++) {
        const GeneId globalGeneId = clusterGraph.geneSet[j];
        const double aExpression = setsAverageExpression[0][j];
        const double bExpression = setsAverageExpression[1][j];
        html <<
            "<tr><td>" << geneNames[globalGeneId] <<
            "<td>" << aExpression <<
            "<td>" << bExpression <<
            "<td>" << max(aExpression, bExpression) <<
            "<td>" << aExpression / bExpression <<
            "<td>" << bExpression / aExpression;
            ;
    }
    // html.precision(oldPrecision);
    html.setf(oldOptions);
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#expressionTable').tablesorter();});"
        "</script>"
        ;
}

