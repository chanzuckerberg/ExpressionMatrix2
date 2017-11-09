// This file contains the implementation of http server functionality
// of class ExpressionMatrix related to cell similarity graphs.

#include "ExpressionMatrix.hpp"
#include "CellGraph.hpp"
#include "color.hpp"
#include "SimilarPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/graph/iteration_macros.hpp>



void ExpressionMatrix::exploreCellGraphs(
    const vector<string>& request,
    ostream& html)
{
    html << "<h1>Cell graphs</h1>";



    // Table with existing cell graphs.
    html <<
        "<table><tr>"
        "<th class=centered>Celll<br>graph<br>name"
        "<th class=centered>Cell<br>set<br>name"
        "<th class=centered>Similar<br>pairs<br>name"
        "<th class=centered>Similarity<br>threshold"
        "<th class=centered>Maximum<br>connectivity"
        "<th class=centered>Number<br>of<br>vertices"
        "<th class=centered>Number<br>of<br>edges"
        "<th class=centered>Number<br>of<br>isolated<br>vertices<br>removed"
        "<th class=centered>Action";
    for(const auto& p: cellGraphs) {
        const string& graphName = p.first;
        const CellGraphInformation& info = p.second.first;
        // const CellGraph& graph = *(p.second.second);
        html << "<tr><td><a href='cellGraph?graphName=" << urlEncode(graphName);
        if(info.edgeCount>50000) {
            // The graph has lots of edges. Don't display them initially.
            html << "&hideEdges=on";
        }
        html << "'>" << graphName << "</a>";
        html << "<td><a href='cellSet?cellSetName=" << urlEncode(info.cellSetName);
        html << "'>" << info.cellSetName << "</a>";
        html << "<td>" << info.similarPairsName;
        html << "<td class=centered>" << info.similarityThreshold;
        html << "<td class=centered>" << info.maxConnectivity;
        html << "<td class=centered>" << info.vertexCount;
        html << "<td class=centered>" << info.edgeCount;
        html << "<td class=centered>" << info.isolatedVertexCount;
        html << "<td class=centered><form action=removeCellGraph><input type=text hidden name=graphName value='" << graphName << "'><input type=submit value='Remove graph " << graphName << "'></form>";
    }



    // Add a final row to the table, containing a form to create a new graph.
    vector<string> availableSimilarPairs;
    getAvailableSimilarPairs(availableSimilarPairs);

    html <<
        "<tr>"
        "<form action=createCellGraph>"
        "<td class=centered><input type=text size=8 required autofocus name=graphName>"
        "<td class=centered>";

    writeCellSetSelection(html, "cellSetName", false);

    html <<
        "<td class=centered><select name=similarPairsName>";
    for(const string& s: availableSimilarPairs) {
        html << "<option value=" << s << ">" << s << "</option>";
    }
    html << "</select>";

    html <<
        "<td class=centered><input type=text style='text-align:center' size=8 name=similarityThreshold value='0.5'>"
        "<td class=centered><input type=text style='text-align:center' size=8 name=maxConnectivity value='20'>"
        "<td><td><td><td class=centered><input type=submit value='Create a new cell graph'>"
        "</form>";

    html << "</table>";



    // Form to compare two graphs.
    if(cellGraphs.size() > 1) {
        html << "<p><form action=compareCellGraphs><input type=submit value=Compare> graphs ";
        writeCellGraphSelection(html, "graphName0", false);
        html << " and ";
        writeCellGraphSelection(html, "graphName1", false);
        html << ".</form>";
    }




}



void ExpressionMatrix::compareCellGraphs(
    const vector<string>& request,
    ostream& html)
{
    // Get the names of the two graphs to be compared.
    string graphName0;
    getParameterValue(request, "graphName0", graphName0);
    string graphName1;
    getParameterValue(request, "graphName1", graphName1);

    // Locate the graphs.
    const auto it0 = cellGraphs.find(graphName0);
    const auto it1 = cellGraphs.find(graphName1);
    if(it0==cellGraphs.end() || it1==cellGraphs.end()) {
        html << "<p>Did not find one or both of the cell graphs to be compared.";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }
    if(it0 == it1) {
        html << "<p>Comparison of a cell graph with itself requested. ";
        html << "<p><form action=cellGraphs><input type=submit value=Continue></form>";
        return;
    }
    const CellGraphInformation& graphCreationParameters0 = it0->second.first;
    const CellGraphInformation& graphCreationParameters1 = it1->second.first;
    const auto graphPointer0 = it0->second.second;
    const auto graphPointer1 = it1->second.second;
    CZI_ASSERT(graphPointer0);
    CZI_ASSERT(graphPointer1);
    const CellGraph& graph0 = *graphPointer0;
    const CellGraph& graph1 = *graphPointer1;



    // Find the common vertices (vertices that correspond to the same cell).
    vector<CellId> cells0, cells1;
    BGL_FORALL_VERTICES(v, graph0, CellGraph) {
        cells0.push_back(graph0[v].cellId);
    }
    BGL_FORALL_VERTICES(v, graph1, CellGraph) {
        cells1.push_back(graph1[v].cellId);
    }
    sort(cells0.begin(), cells0.end());
    sort(cells1.begin(), cells1.end());
    vector<CellId> commonCells;
    std::set_intersection(cells0.begin(), cells0.end(), cells1.begin(), cells1.end(), back_inserter(commonCells));



    // Maps of the edges. Keyed by pair(CellId, CellId), with the lowest numbered cell first.
    // Values: similarities.
    map< pair<CellId, CellId>, float> edgeMap0, edgeMap1;
    BGL_FORALL_EDGES(e, graph0, CellGraph) {
        const CellGraph::vertex_descriptor vA = source(e, graph0);
        const CellGraph::vertex_descriptor vB = target(e, graph0);
        CellId cellIdA = graph0[vA].cellId;
        CellId cellIdB = graph0[vB].cellId;
        CZI_ASSERT(cellIdA != cellIdB);
        if(cellIdB < cellIdA) {
            swap(cellIdA, cellIdB);
        }
        CZI_ASSERT(cellIdA < cellIdB);
        edgeMap0.insert(make_pair( make_pair(cellIdA, cellIdB), graph0[e].similarity));
    }
    BGL_FORALL_EDGES(e, graph1, CellGraph) {
        const CellGraph::vertex_descriptor vA = source(e, graph1);
        const CellGraph::vertex_descriptor vB = target(e, graph1);
        CellId cellIdA = graph1[vA].cellId;
        CellId cellIdB = graph1[vB].cellId;
        CZI_ASSERT(cellIdA != cellIdB);
        if(cellIdB < cellIdA) {
            swap(cellIdA, cellIdB);
        }
        CZI_ASSERT(cellIdA < cellIdB);
        edgeMap1.insert(make_pair( make_pair(cellIdA, cellIdB), graph1[e].similarity));
    }



    // Count edges in common between the two graphs (same cells), and edges
    // that in addition to having the same cells have the same exact similarity.
    // Also, compute the distribution of the discrepancy of the edge similarities
    // for common edges.
    size_t commonEdgeCount = 0;
    size_t commonIdenticalEdgeCount = 0;
    map<float, size_t> similarityDeltaHistogram;    // Key: rounded similarity discrepancy
    const float binWidth = 0.01f;
    for(const auto& p0: edgeMap0) {
        const auto it1 = edgeMap1.find(p0.first);
        if(it1 != edgeMap1.end()) {
            ++commonEdgeCount;
            const float& similarity0 = p0.second;
            const float& similarity1 = it1->second;
            if(similarity0 == similarity1) {
                ++commonIdenticalEdgeCount;
            }
            const float delta = similarity1 - similarity0;
            const float roundedDelta = binWidth * std::round(delta / binWidth);
            const auto it = similarityDeltaHistogram.find(roundedDelta);
            if(it == similarityDeltaHistogram.end()) {
                similarityDeltaHistogram.insert(make_pair(roundedDelta, 1));
            } else {
                ++(it->second);
            }
        }
    }



    // Write a header.
    html << "<h1>Compare cell similarity graphs " << graphName0 << " and " << graphName1 << "</h1>";

    // Begin the table containing the graph comparison.
    // The last column contains checkmarks when the two graphs are in agreement.
    html << "<table><tr><th><th>" << graphName0 << "<th>" << graphName1 << "<th style='width:30px;'>";
    const string checkMark = "&#x2714"; // 2713=checkmark, 2714=heavy checkmark

    // Cell sets.
    html << "<tr><td>Cell set name";
    html << "<td><a href='cellSet?cellSetName=" << urlEncode(graphCreationParameters0.cellSetName);
    html << "'>" << graphCreationParameters0.cellSetName << "</a>";
    html << "<td><a href='cellSet?cellSetName=" << urlEncode(graphCreationParameters1.cellSetName);
    html << "'>" << graphCreationParameters1.cellSetName << "</a>";
    html << "<td class=centered>";
    if(graphCreationParameters0.cellSetName == graphCreationParameters1.cellSetName) {
        html << checkMark;
    }


    // Similar pairs.
    html << "<tr><td>Similar pairs name";
    html << "<td>" << graphCreationParameters0.similarPairsName;
    html << "<td>" << graphCreationParameters1.similarPairsName;
    html << "<td class=centered>";
    if(graphCreationParameters0.similarPairsName == graphCreationParameters1.similarPairsName) {
        html << checkMark;
    }

    // Similarity thresholds.
    html << "<tr><td>Similarity threshold";
    html << "<td class=centered>" << graphCreationParameters0.similarityThreshold;
    html << "<td class=centered>" << graphCreationParameters1.similarityThreshold;
    html << "<td class=centered>";
    if(graphCreationParameters0.similarityThreshold == graphCreationParameters1.similarityThreshold) {
        html << checkMark;
    }

    // Maximum connectivities.
    html << "<tr><td>Maximum connectivity";
    html << "<td class=centered>" << graphCreationParameters0.maxConnectivity;
    html << "<td class=centered>" << graphCreationParameters1.maxConnectivity;
    html << "<td class=centered>";
    if(graphCreationParameters0.maxConnectivity == graphCreationParameters1.maxConnectivity) {
        html << checkMark;
    }

    // Number of vertices.
    const size_t vertexCount0 = cells0.size();
    const size_t vertexCount1 = cells1.size();
    html << "<tr><td>Number of vertices (cells)";
    html << "<td class=centered>" << vertexCount0;
    html << "<td class=centered>" << vertexCount1;
    html << "<td class=centered>";
    if(vertexCount0 == vertexCount1) {
        html << checkMark;
    }

    // Number of isolated vertices.
    const size_t isolatedVertexCount0 = graphCreationParameters0.isolatedVertexCount;
    const size_t isolatedVertexCount1 = graphCreationParameters1.isolatedVertexCount;
    html << "<tr><td>Number of isolated vertices (cells);";
    html << "<td class=centered>" << isolatedVertexCount0;
    html << "<td class=centered>" << isolatedVertexCount1;
    html << "<td class=centered>";
    if(isolatedVertexCount0 == isolatedVertexCount1) {
        html << checkMark;
    }

    // Number of vertices in common between the two graphs.
    html << "<tr><td>Number of common vertices (cells)";
    html << "<td class=centered>" << commonCells.size();
    html << "<td class=centered>" << commonCells.size();
    html << "<td class=centered>";
    if(vertexCount0==commonCells.size() && vertexCount1==commonCells.size()) {
        html << checkMark;
    }

    // Number of vertices not in common between the two graphs.
    html << "<tr><td>Number of vertices (cells) not in common";
    html << "<td class=centered>" << vertexCount0-commonCells.size();
    html << "<td class=centered>" << vertexCount1-commonCells.size();
    html << "<td class=centered>";
    if(vertexCount0==commonCells.size() && vertexCount1==commonCells.size()) {
        html << checkMark;
    }

    // Number of edges.
    const size_t edgeCount0 = edgeMap0.size();
    const size_t edgeCount1 = edgeMap1.size();
    html << "<tr><td>Number of edges";
    html << "<td class=centered>" << edgeCount0;
    html << "<td class=centered>" << edgeCount1;
    html << "<td class=centered>";
    if(edgeCount0 == edgeCount1) {
        html << checkMark;
    }

    // Number of common edges (same cells).
    html << "<tr><td>Number of common edges (same cells)";
    html << "<td class=centered>" << commonEdgeCount;
    html << "<td class=centered>" << commonEdgeCount;
    html << "<td class=centered>";
    if(edgeCount0==commonEdgeCount && edgeCount1==commonEdgeCount) {
        html << checkMark;
    }

    // Number of identical common edges (same cells and same similarity).
    html << "<tr><td>Number of identical common edges (same cells and same similarity)";
    html << "<td class=centered>" << commonIdenticalEdgeCount;
    html << "<td class=centered>" << commonIdenticalEdgeCount;
    html << "<td class=centered>";
    if(edgeCount0==commonIdenticalEdgeCount && edgeCount1==commonIdenticalEdgeCount) {
        html << checkMark;
    }

    // End the table containing the graph comparison.
    html << "</table>";



    // Table with the distribution of the discrepancy in similarity for the common edges.
    html <<
        "<h2>Distribution of similarity difference</h2>"
        "<p>This table shows the distribution of the difference between the similarity "
        " for an edge of graph " << graphName1 << " and the similarity for the corresponding edge of graph " << graphName0 <<
        " computed over all " << commonEdgeCount << " common edges (same cells)."
        "<table><tr><th>Difference<th>Frequency";
    for(const auto& p: similarityDeltaHistogram) {
        html << "<tr><td class=centered>" << p.first << "<td class=centered>" << p.second;
    }
    html << "</table>";
}




void ExpressionMatrix::exploreCellGraph(
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

    // Find the graph.
    const auto it = cellGraphs.find(graphName);
    if(it == cellGraphs.end()) {
        html << "<p>Graph " << graphName << " does not exists.";
        return;
    }
    const CellGraphInformation& graphInformation = it->second.first;
    const string& similarPairsName = graphInformation.similarPairsName;
    vector<string> geneSetNames = geneSetNamesFromSimilarPairsName(similarPairsName);
    const string geneSetName = geneSetNames.empty() ? "" : geneSetNames.front();
    CellGraph& graph = *(it->second.second);

    // Write the title.
    html << "<h1>Graph " << graphName << "</h1>";


    // Compute the graph layout, if necessary.
    if (!graph.layoutWasComputed) {
        html << "<div style='font-family:courier'>";
        html << timestamp << "Graph layout computation begins.";
        graph.computeLayout();
        html << "<br>" << timestamp << "Graph layout computation ends.";
        html << "</div>";
        graph.layoutWasComputed = true;
    }


    // Div to contain the table and the form for coloring.
    html << "<div>";

    // Write a table with the graph creation parameters and other information.
    html << "<div style='float:left;margin:10px'>";
    html << "<table>";
    html << "<tr><td>Cell set name<td><a href='cellSet?cellSetName=" << urlEncode(graphInformation.cellSetName);
    html << "'>" << graphInformation.cellSetName << "</a>";
    html << "<tr><td>Similar pairs name<td>" << graphInformation.similarPairsName;
    html << "<tr><td>Gene set name<td>" << geneSetName;
    html << "<tr><td>Similarity threshold<td class=centered>" << graphInformation.similarityThreshold;
    html << "<tr><td>Maximum connectivity<td class=centered>" << graphInformation.maxConnectivity;
    html << "<tr><td>Number of vertices (cells)<td class=centered>" << boost::num_vertices(graph);
    html << "<tr><td>Number of edges<td class=centered>" << boost::num_edges(graph);
    html << "<tr><td>Number of isolated vertices (cells) removed<td class=centered>" << graphInformation.isolatedVertexCount;
    html << "</table>";
    html << "</div>";




    // Compute maximum and minimum coordinates for the vertices in the graph,
    // then enlarge them to a square.
    double xMin, xMax, yMin, yMax;
    graph.computeCoordinateRange(xMin, xMax, yMin, yMax);
    const double deltaX = xMax - xMin;
    const double deltaY = yMax - yMin;
    double delta;
    if(deltaX <= deltaY) {
        delta = deltaY;
        xMax = xMin + delta;
    } else {
        delta = deltaX;
        yMax = yMin + delta;
    }
    delta *= 1.05;  // Leave a bit of space around it.

    // Extract parameters from the request.
    double svgSizePixels = 600;
    getParameterValue(request, "svgSizePixels", svgSizePixels);
    double vertexRadius = 1.5 * delta / svgSizePixels;
    getParameterValue(request, "vertexRadius", vertexRadius);
    double edgeThickness = 0.05 * delta / svgSizePixels;
    getParameterValue(request, "edgeThickness", edgeThickness);
    double xViewBoxCenter = (xMin + xMax) / 2.;
    getParameterValue(request, "xViewBoxCenter", xViewBoxCenter);
    double yViewBoxCenter = (yMin + yMax) / 2.;
    getParameterValue(request, "yViewBoxCenter", yViewBoxCenter);
    double viewBoxHalfSize = delta / 2.;
    getParameterValue(request, "viewBoxHalfSize", viewBoxHalfSize);
    string geneIdStringForColoringByGeneExpression;
    getParameterValue(request, "geneIdForColoringByGeneExpression", geneIdStringForColoringByGeneExpression);
    const NormalizationMethod normalizationMethod = getNormalizationMethod(request, NormalizationMethod::L2);
    string cellIdStringForColoringBySimilarity;
    getParameterValue(request, "cellIdStringForColoringBySimilarity", cellIdStringForColoringBySimilarity);
    CellId cellIdForColoringBySimilarity = cellIdFromString(cellIdStringForColoringBySimilarity);
    string metaDataName;
    getParameterValue(request, "metaDataName", metaDataName);
    string metaDataMeaning = "category";
    getParameterValue(request, "metaDataMeaning", metaDataMeaning);
    string coloringOption = "noColoring";
    getParameterValue(request, "coloringOption", coloringOption);
    string reuseColors = "off";
    getParameterValue(request, "reuseColors", reuseColors);
    string hideEdges = "off";
    getParameterValue(request, "hideEdges", hideEdges);

    // Values corresponding to the minimum and maximum color.
    double minColorValue, maxColorValue;
    const double minColorValueIsPresent = getParameterValue(request, "minColorValue", minColorValue);
    const double maxColorValueIsPresent = getParameterValue(request, "maxColorValue", maxColorValue);



    // Form used to specify graph display options.
    // It has an onsubmit function that blanks out minColorValue and maxColorValue
    // before submission if any options affecting the validity of minColorValue and maxColorValue
    // have changed.
    html <<
        "<script>"
        "function massageForm() {\n"
        "    coloringOptionChanged = !document.getElementById('" << coloringOption << "Radio').checked;\n"
        "    // window.alert('Coloring option changed: ' + coloringOptionChanged);\n"
        "    geneIdChanged = document.getElementById('geneIdInput').value != '" << geneIdStringForColoringByGeneExpression << "';\n"
        "    // window.alert('Gene id changed: ' + geneIdChanged);\n"
        "    normalizationMethodChanged = document.getElementById('normalizationMethod').value != '" <<
        normalizationMethodToShortString(normalizationMethod) << "';\n"
        "    // window.alert('Normalization method changed: ' + normalizationMethodChanged);\n"
        "    cellIdChanged = document.getElementById('cellIdInput').value != '" << cellIdStringForColoringBySimilarity << "';\n"
        "    // window.alert('Cell id changed: ' + cellIdChanged);\n"
        "    metaDataChanged = document.getElementById('metaDataName').value != '" << metaDataName << "';\n"
        "    // window.alert('Meta data changed: ' + metaDataChanged);\n"
        "    if(coloringOptionChanged || geneIdChanged || normalizationMethodChanged || cellIdChanged || metaDataChanged) {\n"
        "        document.getElementById('minColorInput').value = '';"
        "        document.getElementById('maxColorInput').value = '';"
        "    }\n"
        "}\n"
        "</script>"
        ;
    html << "<div><form id=coloringForm onsubmit='massageForm()'>";



    // Radio button to request no coloring.
    html <<
        "<br><input type=radio name=coloringOption value=noColoring id=noColoringRadio";
    if(coloringOption == "noColoring" ||
        (coloringOption!="bySimilarity" && coloringOption!="byCluster" && coloringOption!="byMetaData")) {
        html << " checked=checked";
    }
    html <<">No coloring";



    // Radio button to color by expression of a specified gene and related options.
    html <<
        "<br><input type=radio name=coloringOption value=byGeneExpression id=byGeneExpressionRadio";
    if(coloringOption == "byGeneExpression") {
        html << " checked=checked";
    }
    html <<
        ">Color by expression of gene "
        "<input type=text name=geneIdForColoringByGeneExpression id=geneIdInput size=8";
    if(!geneIdStringForColoringByGeneExpression.empty()) {
        html << " value='" << geneIdStringForColoringByGeneExpression << "'";
    }
    html << ">";

    html << " using ";
    writeNormalizationSelection(html, normalizationMethod);




    // Radio button to color by similarity to a specified cell and related options.
    html <<
        "<br><input type=radio name=coloringOption value=bySimilarity id=bySimilarityRadio";
    if(coloringOption == "bySimilarity") {
        html << " checked=checked";
    }
    html <<
        ">Color by similarity to cell "
        "<input type=text name=cellIdStringForColoringBySimilarity id=cellIdInput size=8";
    if(!cellIdStringForColoringBySimilarity.empty()) {
        html << " value='" << cellIdStringForColoringBySimilarity << "'";
    }
    html << ">";



    // Radio button to color by meta data and related options.
    html <<
        "<br><input type=radio name=coloringOption value=byMetaData id=byMetaDataRadio";
    if(coloringOption == "byMetaData") {
        html << " checked=checked";
    }
    html <<">Color by meta data field ";
    set<string> selectedMetaData;
    if(!metaDataName.empty()) {
        selectedMetaData.insert(metaDataName);
    }
    writeMetaDataSelection(html, "metaDataName", selectedMetaData, false);
    html << " interpreted as a <select name=metaDataMeaning>";
    html << "<option value=category";
    if(metaDataMeaning == "category") {
        html << " selected=selected";
    }
    html << ">category</option>";
    html << "<option value=number";
    if(metaDataMeaning == "number") {
        html << " selected=selected";
    }
    html << ">number</option>";
    html << "<option value=color";
    if(metaDataMeaning == "color") {
        html << " selected=selected";
    }
    html << ">color</option>";
    html << "</select>";
    html << ". Allow using a color for multiple non-adjacent categories <input type=checkbox name=reuseColors";
    if(reuseColors=="on") {
        html << " checked";
    }
    html << ">";
    html << "<br>Hide graph edges  <input type=checkbox name=hideEdges";
    if(hideEdges=="on") {
        html << " checked";
    }
    html << ">";



    // Add a few necessary hidden input fields, so geometric parameters maintained in Javascript code
    // can be remembered when the user clicks "Redraw graph.".
    html <<
        "<input type=hidden name=graphName value='" << graphName << "'>" <<
        "<input id=svgSizePixelsInput type=hidden name=svgSizePixels value='" << svgSizePixels << "'>"
        "<input id=vertexRadiusInput type=hidden name=vertexRadius value='" << vertexRadius << "'>"
        "<input id=edgeThicknessInput type=hidden name=edgeThickness value='" << edgeThickness << "'>"
        "<input id=xViewBoxCenterInput type=hidden name=xViewBoxCenter value='" << xViewBoxCenter << "'>"
        "<input id=yViewBoxCenterInput type=hidden name=yViewBoxCenter value='" << yViewBoxCenter << "'>"
        "<input id=viewBoxHalfSizeInput type=hidden name=viewBoxHalfSize value='" << viewBoxHalfSize << "'>"
        ;



    // Add a submit button and finish the form.
    html <<
        "<br><button type=submit>Redraw graph</button>"
        "</form></div>";


    // End the div to contain the table and the form for coloring.
    html << "</div>";


    // Alternatives for the arrow symbols used below.
    // right        left        up          down
    // &#129094;    &#129092;   &#129093;   &#129095;   // Original fat arrows, don't work on Mac Chrome
    // &#10145;     &#11013;    &#11014;    &#11015;    // Suggested by Andrey, less fat but seem to work everywhere
    // &#8594;      &#8592;     &#8593;     &#8595;     // Suggested by W3C, but they are not fat
    // &rarr;       &larr;      &uarr;      &darr;      // Same as above, using HTML symbols


    // Controls and code to change the graph display.
    html <<
        "<script>"
        "var xViewBoxCenter=" << xViewBoxCenter << ";"
        "var yViewBoxCenter=" << yViewBoxCenter << ";"
        "var viewBoxHalfSize=" << viewBoxHalfSize << ";"
        "var scaleFactor = 1.1;"
        "var largerFactor = 1.1;"
        "var translationDelta = 0.2;"
        "var widthAndHeight = " << svgSizePixels << ";"
        "var vertexRadius = " << vertexRadius << ";"
        "var vertexRadiusFactor = 1.2;"
        "var edgeThickness = " << edgeThickness << ";"
        "var edgeThicknessFactor = 1.2;"
        ;
    html << R"%(
    function rescaleSvg()
    {
        var svgElement = document.getElementById("graphSvg");
        var xMin = xViewBoxCenter - viewBoxHalfSize;
        var yMin = yViewBoxCenter - viewBoxHalfSize;
        var width = 2 * viewBoxHalfSize;
        var height = 2 * viewBoxHalfSize;
        var viewBox = xMin + " " + yMin + " " + width + " " + height;
        svgElement.setAttribute("viewBox", viewBox);
        document.getElementById("xViewBoxCenterInput").value = xViewBoxCenter;
        document.getElementById("yViewBoxCenterInput").value = yViewBoxCenter;
        document.getElementById("viewBoxHalfSizeInput").value = viewBoxHalfSize;
    }
    function zoomIn()
    {
        viewBoxHalfSize /= scaleFactor;
        rescaleSvg();
    }
    function zoomOut()
    {
        viewBoxHalfSize *= scaleFactor;
        rescaleSvg();
    }
    function moveRight()
    {
        xViewBoxCenter += viewBoxHalfSize * translationDelta;
        rescaleSvg();
    }
    function moveLeft()
    {
        xViewBoxCenter -= viewBoxHalfSize * translationDelta;
        rescaleSvg();
    }
    function moveUp()
    {
        yViewBoxCenter -= viewBoxHalfSize * translationDelta;
        rescaleSvg();
    }
    function moveDown()
    {
        yViewBoxCenter += viewBoxHalfSize * translationDelta;
        rescaleSvg();
    }
    function resizeSvg()
    {
        var svgElement = document.getElementById("graphSvg");
        widthAndHeight = Math.round(widthAndHeight);
        svgElement.setAttribute("width", widthAndHeight);
        svgElement.setAttribute("height", widthAndHeight);
        document.getElementById("svgSizePixelsInput").value = widthAndHeight;
    }
    function larger()
    {
        widthAndHeight *= largerFactor;
        resizeSvg();
    }
    function smaller()
    {
        widthAndHeight /= largerFactor;
        resizeSvg();
    }
    function updateVertexRadius()
    {
        var vertexGroups = document.getElementById("vertices").childNodes; 
        for(i=0; i<vertexGroups.length; i++) {
             var groupVertices = vertexGroups[i].childNodes;
             for(j=0; j<groupVertices.length; j++) {
                 // Each vertex is inside an <a> element, so we have to do childNodes[0] to access it.
                 groupVertices[j].childNodes[0].setAttribute("r", vertexRadius);
             }
        }
        document.getElementById("vertexRadiusInput").value = vertexRadius;
    }
    function increaseVertexRadius()
    {
        vertexRadius *= vertexRadiusFactor;
        updateVertexRadius();
    }
    function decreaseVertexRadius()
    {
        vertexRadius /= vertexRadiusFactor;
        updateVertexRadius();
    }
    function updateEdgeThickness()
    {
        var edges = document.getElementById("edges").childNodes; 
        for(i=0; i<edges.length; i++) {
             edges[i].style.strokeWidth = edgeThickness;
        }
        document.getElementById("edgeThicknessInput").value = edgeThickness;
    }
    function increaseEdgeThickness()
    {
        edgeThickness *= edgeThicknessFactor;
        updateEdgeThickness();
    }
    function decreaseEdgeThickness()
    {
        edgeThickness /= edgeThicknessFactor;
        updateEdgeThickness();
    }
</script>
<div style='clear:both;'>
<button onClick='larger()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
Larger
</button>
<button onClick='smaller()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
Smaller
</button>
<button onClick='zoomIn()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
Zoom in 
</button>
<button onClick='zoomOut()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
Zoom out 
</button>
<br>
<button onClick='increaseVertexRadius()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
Larger vertex
</button>
<button onClick='decreaseVertexRadius()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
Smaller vertex
</button>
<button onClick='increaseEdgeThickness()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
Thicker edge
</button>
<button onClick='decreaseEdgeThickness()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
Thinner edge
</button>
<br>
<button onClick='moveRight()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
&#10145;
</button>
<button onClick='moveLeft()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
&#11013;
</button>
<button onClick='moveUp()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
&#11014;
</button>
<button onClick='moveDown()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
&#11015;
</button>
</div>
    )%";


    // Some data structures that need to be defined at this level.
    map<string, int> groupMap;                          // Maps meta data string to group number.
    map<int, string> colorMap;                          // Maps group number to color string.
    vector< pair<int, string> > sortedFrequencyTable;   // Pairs (frequency, meta data string)
    int couldNotColor = 0;

    // Flag that will be set to true if we are coloring by number, that is, using a continuous scale.
    // In the case we store in each vertex the value that we want to color by.
    bool colorByNumber = false;



    // Color the graph by expression of a given gene.
    if(coloringOption == "byGeneExpression") {
        const GeneId geneId = geneIdFromString(geneIdStringForColoringByGeneExpression);
        if(geneId == invalidGeneId) {
            html << "<p>Gene not found.";
            return;
        }
        colorByNumber = true;

        // Set the value field for all the vertices.
        BGL_FORALL_VERTICES(v, graph, CellGraph) {
            CellGraphVertex& vertex = graph[v];
            const double rawCount = getCellExpressionCount(vertex.cellId, geneId);
            if(normalizationMethod == NormalizationMethod::None) {
                vertex.value = rawCount;
            } else {
                const Cell& cell = cells[vertex.cellId];
                if(normalizationMethod == NormalizationMethod::L1) {
                    vertex.value = rawCount * cell.norm1Inverse;
                } else if(normalizationMethod == NormalizationMethod::L2) {
                    vertex.value = rawCount * cell.norm2Inverse;
                } else if(normalizationMethod == NormalizationMethod::Invalid){
                    html << "<p>Invalid normalization method.";
                    return;
                }
            }
        }
    }


    // Color the graph by similarity to a specified cell.
    else if(coloringOption == "bySimilarity") {
        if(cellIdForColoringBySimilarity<0 || size_t(cellIdForColoringBySimilarity) >= cells.size()) {
            html << "<p>Invalid cell id.";
            return;
        }
        colorByNumber = true;
        const SimilarPairs similarPairs(directoryName + "/SimilarPairs-" + similarPairsName, true);
        const GeneSet& geneSet = similarPairs.getGeneSet();
        BGL_FORALL_VERTICES(v, graph, CellGraph) {
            CellGraphVertex& vertex = graph[v];
            vertex.value = computeCellSimilarity(geneSet, cellIdForColoringBySimilarity, vertex.cellId);
        }
    }



    // Color the graph by metadata.
    // Each vertex receives a color determined by the chosen meta data field.
    // If interpretMetaDataAsColor is "on", the meta data is interpreted directly as an html color.
    // Otherwise, it is interpreted as a category and mapped to a color.
    else if(coloringOption == "byMetaData") {



        // Color the graph by meta data, interpreting the meta data as a category.
        if(metaDataMeaning == "category") {

            // We need to assign groups based on the of values of the specified meta data field.
            // Find the frequency of each of them.
            map<string, int> frequencyTable;
            BGL_FORALL_VERTICES(v, graph, CellGraph) {
                const CellGraphVertex& vertex = graph[v];
                const string metaDataValue = getCellMetaData(vertex.cellId, metaDataName);
                const auto it = frequencyTable.find(metaDataValue);
                if(it == frequencyTable.end()) {
                    frequencyTable.insert(make_pair(metaDataValue, 1));
                } else {
                    ++(it->second);
                }
            }

            // Sort them by decreasing frequency.
            for(const auto& p: frequencyTable) {
                sortedFrequencyTable.push_back(make_pair(p.second, p.first));
            }
            sort(sortedFrequencyTable.begin(), sortedFrequencyTable.end(), std::greater< pair<int, string> >());

            // Map the meta data categories to groups.
            for(size_t group=0; group<sortedFrequencyTable.size(); group++) {
                groupMap.insert(make_pair(sortedFrequencyTable[group].second, group));
            }

            // Assign the vertices to groups..
            BGL_FORALL_VERTICES(v, graph, CellGraph) {
                CellGraphVertex& vertex = graph[v];
                const string metaData = getCellMetaData(vertex.cellId, metaDataName);
                vertex.group = groupMap[metaData];
            }



            // Map the groups to colors.
            if(reuseColors == "on") {

                // Each color can be used for more than one category (group),
                // as long as the vertices of every edge have distinct colors.
                vector<uint32_t> graphColoringTable;
                graph.assignColorsToGroups(graphColoringTable);
                for(size_t group=0; group<sortedFrequencyTable.size(); group++) {
                    const uint32_t iColor = graphColoringTable[group];
                    string colorString = "black";
                    if(iColor < 12) {
                        colorString = colorPalette1(iColor);
                    }
                    colorMap.insert(make_pair(group, colorString));
                }

            } else {

                // Each color gets used for a single meta data category.
                for(size_t group=0; group<sortedFrequencyTable.size(); group++) {
                    string colorString = "black";
                    if(group<12) {
                        colorString = colorPalette1(group);
                    } else {
                        couldNotColor += sortedFrequencyTable[group].first;
                    }
                    colorMap.insert(make_pair(group, colorString));
                }
            }

        }



        // Color the graph by meta data, interpreting the meta data as an html color
        //  (that is, color name, or # followed by 6 hex digits).
        else if(metaDataMeaning == "color") {

            // The meta data field is interpreted directly as an html color.
            BGL_FORALL_VERTICES(v, graph, CellGraph) {
                CellGraphVertex& vertex = graph[v];
                vertex.group = 0;
                vertex.color = getCellMetaData(vertex.cellId, metaDataName);
            }
        }



        // Color by meta data, interpreting the meta data value as a number.
        // We store in each vertex the meta data value that will determine the vertex color.
        else if(metaDataMeaning == "number") {
            colorByNumber = true;
            BGL_FORALL_VERTICES(v, graph, CellGraph) {
                CellGraphVertex& vertex = graph[v];
                vertex.group = 0;
                vertex.value = std::numeric_limits<double>::max();
                try {
                    vertex.value = lexical_cast<double>(getCellMetaData(vertex.cellId, metaDataName));
                } catch(bad_lexical_cast) {
                    // If the meta data cannot be interpreted as a number, the value is left at
                    // the ":invalid" value set above, and the vertex will be colored black.
                }
                // The vertex color will be computed later, so the code can be shared with other
                // coloring options.
            }
        }


        // Otherwise, don't color the vertices.
        else {
            BGL_FORALL_VERTICES(v, graph, CellGraph) {
                CellGraphVertex& vertex = graph[v];
                vertex.color.clear();
                vertex.group = 0;
            }
        }



        // Don't color the edges.
        BGL_FORALL_EDGES(e, graph, CellGraph) {
            graph[e].color.clear();
        }

    }



    // Otherwise, all vertices and edges are colored black.
    else {

        BGL_FORALL_VERTICES(v, graph, CellGraph) {
            CellGraphVertex& vertex = graph[v];
            vertex.color.clear();
            vertex.group = 0;
        }
        BGL_FORALL_EDGES(e, graph, CellGraph) {
            graph[e].color.clear();
        }
    }



    // If coloring by number, compute the color of each vertex.
    double minValue = std::numeric_limits<double>::max();
    double maxValue = std::numeric_limits<double>::lowest();
    if(colorByNumber) {

        // Compute the minimum and maximum values.
        BGL_FORALL_VERTICES(v, graph, CellGraph) {
            const double value = graph[v].value;
            if(value == std::numeric_limits<double>::max()) {
                continue;
            }
            minValue = min(minValue, value);
            maxValue = max(maxValue, value);
        }

        if(!minColorValueIsPresent) {
            minColorValue = minValue;
        }
        if(!maxColorValueIsPresent) {
            maxColorValue = maxValue;
        }

        // Now compute the colors.
        if(minValue==maxValue || maxValue==std::numeric_limits<double>::lowest()) {
            BGL_FORALL_VERTICES(v, graph, CellGraph) {
                graph[v].color = "black";
            }
        } else {
            const double scalingFactor = 1./(maxColorValue - minColorValue);
            BGL_FORALL_VERTICES(v, graph, CellGraph) {
                CellGraphVertex& vertex = graph[v];
                const double value = vertex.value;
                if(value == std::numeric_limits<double>::max()) {
                    continue;
                }
                vertex.color = spectralColor(scalingFactor * (value-minColorValue));
            }
        }
    }



    // Begin a div to contain the graphics and the table of meta data groups
    // (or other information depending on the coloring option used), side by side.
    html << "<div>";

    // Write the graph as an svg object.
    html << "<div style='float:left;margin:10px''>";
    graph.writeSvg(html, hideEdges=="on", svgSizePixels, xViewBoxCenter, yViewBoxCenter, viewBoxHalfSize, vertexRadius, edgeThickness, colorMap, geneSetName);
    html << "</div>";



    // Additional processing for coloring by meta data interpreted as category.
    if(coloringOption=="byMetaData" && metaDataMeaning=="category") {


        // Write the table of meta data groups.
        html << "<div><table>";
        if(reuseColors != "on") {
            html << "<tr><th>" << metaDataName << "<th>Frequency<th>Color";
            for(size_t color=0; color<min(size_t(12), sortedFrequencyTable.size()); color++) {
                const auto& p = sortedFrequencyTable[color];
                html << "<tr id=vertexGroupRow" << color << "><td>" << p.second << "<td class=centered>" << p.first <<
                    "<td  class=centered style='width:20px;background-color:" << colorMap[groupMap[p.second]] << "'>";
            }
            if(sortedFrequencyTable.size() > 12) {
                html << "<tr><td>All others<td class=centered>" << couldNotColor <<
                    "<td  class=centered style='width:20px;background-color:black'>";
            }
        }
        html << "<tr id=highlightedMetaDataRow><td id=highlightedMetaData><td colspan=2>Currently highlighted";
        html << "</table>";

        // Code to highlight groups
        html << R"%(
<script>
function highlight(name)
{
    element = document.getElementById(name);
    element.style.oldFill = element.style.fill;
    element.style.fill = 'red';
}
function unhighlight(name)
{
    element = document.getElementById(name);
    element.style.fill = element.style.oldFill;
}
function writeHighlightedMetaData(name)
{
    document.getElementById('highlightedMetaData').innerHTML = name;
    document.getElementById('highlightedMetaDataRow').style.backgroundColor= 'pink';
}
function removeHighlightedMetaData()
{
    document.getElementById('highlightedMetaData').innerHTML = '';
    document.getElementById('highlightedMetaDataRow').style.backgroundColor= 'white';
}
    )%";
    html <<
        "function highlightTableRow(groupNumber) {\n"
        "    document.getElementById('vertexGroupRow' + groupNumber).style.backgroundColor = 'pink';\n"
        "}\n"
        "function unhighlightTableRow(groupNumber) {\n"
        "    document.getElementById('vertexGroupRow' + groupNumber).style.backgroundColor = 'white';\n"
        "}\n"
        "var element;\n";
        for(const auto& p: groupMap) {
            html <<
                "element = document.getElementById('vertexGroup" << p.second << "');\n"
                "element.onmouseover = function(){highlight('vertexGroup" << p.second << "'); writeHighlightedMetaData('" << p.first << "');";
            if(reuseColors != "on") {
                html << " highlightTableRow('" << p.second << "');";
            }
            html <<
                "};\n"
                "element.onmouseout = function(){unhighlight('vertexGroup" << p.second << "'); removeHighlightedMetaData();";
            if(reuseColors != "on") {
                html << " unhighlightTableRow('" << p.second << "');";
            }
            html << "};\n";
        }
        html << "</script>";

        // End the div containing the table.
        html << "</div>";

    }



    // Additional processing for coloring by number.
    if(colorByNumber) {

        if(!(minValue==maxValue || maxValue==std::numeric_limits<double>::lowest())) {

            // Write the color legend
            html << "<div><table>";
            html << "<tr style='height:30px;'><td>Minimum<td class=centered>" << minValue;
            const int n = 4;
            for(int i=0; i<=n; i++) {
                const double x = i * (1./n);
                html << "<tr style='height:30px';><td style='background-color:" << spectralColor(x) << "'>";
                html << "<td class=centered>";
                if(i==0) {
                    html << "<input type=text style='text-align:center;background-color:LightGrey;' form=coloringForm name=minColorValue id=minColorInput value='" << minColorValue << "'>";
                } else if(i==n) {
                    html << "<input type=text style='text-align:center;background-color:LightGrey;' form=coloringForm name=maxColorValue id=maxColorInput value='" << maxColorValue << "'>";
                } else {
                    html << minColorValue + x*(maxColorValue-minColorValue);;
                }
            }
            html << "<tr style='height:30px;'><td>Maximum<td class=centered>" << maxValue;
            html << "</table>";
            html << "You can change the color scale by editing the cells with a grey background.";
        }
    }



    // End the div to contain the graphics and the table of meta data groups, side by side.
    html << "</div>";
}
