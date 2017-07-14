#include "ExpressionMatrix.hpp"
#include "color.hpp"
#include "CellSimilarityGraph.hpp"
#include "ClusterGraph.hpp"
#include "orderPairs.hpp"
#include "randIndex.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "boost_tuple_tuple.hpp"
#include <boost/filesystem.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

#include "fstream.hpp"



// Function that provides simple http functionality
// to facilitate data exploration and debugging.
// It is passed the string of the GET request,
// already parsed using "?=" as separators.
// The request is guaranteed not to be empty.
void ExpressionMatrix::processRequest(
    const vector<string>& request,
    ostream& html)
{
    // Look up the keyword to find the function that will process this request.
    const string& keyword = request.front();
    const auto it = serverFunctionTable.find(keyword);


    // If we did not find the keyword, see if it is a documentation request.
    if(it == serverFunctionTable.end()) {

    	// See if it is of the form help/xyz.
    	// Here, xyz is not allowed to contain any slashes,
    	// otherwise we would open a big security hole: requests of the form help/../../file
    	// would give access to anywhere in the file system.
    	vector<string> tokens;
        boost::algorithm::split(tokens, keyword, boost::algorithm::is_any_of("/"));
        if(tokens.size()==3 && tokens.front().empty() && tokens[1]=="help") {
        	ifstream file(serverParameters.docDirectory + "/" + tokens[2]);
        	if(file) {
        		html << "\r\n" << file.rdbuf();
        		return;
        	}
        }

    	// This is not a documentation request. Write a message and stop here.
        html << "\r\nUnsupported keyword " << keyword;
        return;
    }



    // We found the keyword. Get the function that processes this keyword.
    const auto function = it->second;



    // Write everything that goes before the html body.
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

    // The processing function is only responsible for writing the html body.
    try {
        (this->*function)(request, html);
    } catch(std::exception& e) {
        html << e.what();
    }

    html << "</body>";
    html << "</html>";
}



void ExpressionMatrix::explore(const ServerParameters& serverParameters)
{
	// Store the server parameters.
	this->serverParameters = serverParameters;

	// Validate the docDirectory in the server parameters.
	// If specified and not valid, set it to an empty string.
	if(this->serverParameters.docDirectory.empty()) {
		cout << "A documentation directory was not specified. Documentation will not be available from the server." << endl;
	} else {
		ifstream file(serverParameters.docDirectory + "/index.html");
		if(!file) {
			cout << "Documentation index file " << serverParameters.docDirectory;
			cout << "/index.html could not be open. Documentation will not be available from the server." << endl;
			this->serverParameters.docDirectory.clear();
		}
	}

	// Invoke teh base class.
	HttpServer::explore(serverParameters.port);
}


void ExpressionMatrix::fillServerFunctionTable()
{
    serverFunctionTable[""]                                 = &ExpressionMatrix::exploreSummary;
    serverFunctionTable["/"]                                = &ExpressionMatrix::exploreSummary;
    serverFunctionTable["/index"]                           = &ExpressionMatrix::exploreSummary;
    serverFunctionTable["/exploreHashTableSummary"]         = &ExpressionMatrix::exploreHashTableSummary;
    serverFunctionTable["/gene"]                            = &ExpressionMatrix::exploreGene;
    serverFunctionTable["/cell"]                            = &ExpressionMatrix::exploreCell;
    serverFunctionTable["/compareTwoCells"]                 = &ExpressionMatrix::compareTwoCells;
    serverFunctionTable["/cellSets"]                        = &ExpressionMatrix::exploreCellSets;
    serverFunctionTable["/cellSet"]                         = &ExpressionMatrix::exploreCellSet;
    serverFunctionTable["/createCellSetUsingMetaData"]      = &ExpressionMatrix::createCellSetUsingMetaData;
    serverFunctionTable["/createCellSetIntersectionOrUnion"]= &ExpressionMatrix::createCellSetIntersectionOrUnion;
    serverFunctionTable["/downsampleCellSet"]				= &ExpressionMatrix::downsampleCellSet;
    serverFunctionTable["/removeCellSet"]                   = &ExpressionMatrix::removeCellSet;
    serverFunctionTable["/graphs"]                          = &ExpressionMatrix::exploreGraphs;
    serverFunctionTable["/compareGraphs"]                   = &ExpressionMatrix::compareGraphs;
    serverFunctionTable["/graph"]                           = &ExpressionMatrix::exploreGraph;
    serverFunctionTable["/clusterDialog"]                   = &ExpressionMatrix::clusterDialog;
    serverFunctionTable["/cluster"]                         = &ExpressionMatrix::cluster;
    serverFunctionTable["/createNewGraph"]                  = &ExpressionMatrix::createNewGraph;
    serverFunctionTable["/removeGraph"]                     = &ExpressionMatrix::removeGraph;
    serverFunctionTable["/metaData"]                        = &ExpressionMatrix::exploreMetaData;
    serverFunctionTable["/metaDataHistogram"]               = &ExpressionMatrix::metaDataHistogram;
    serverFunctionTable["/metaDataContingencyTable"]        = &ExpressionMatrix::metaDataContingencyTable;
    serverFunctionTable["/removeMetaData"]                  = &ExpressionMatrix::removeMetaData;
}



void ExpressionMatrix::writeNavigation(ostream& html)
{
    writeNavigation(html, "Summary", "index");
    writeNavigation(html, "Genes", "gene");
    writeNavigation(html, "Cells", "cell");
    writeNavigation(html, "Compare two cells", "compareTwoCells");
    writeNavigation(html, "Cell sets", "cellSets");
    writeNavigation(html, "Cell meta data", "metaData");
    writeNavigation(html, "Graphs", "graphs");

    if(!serverParameters.docDirectory.empty()) {
        writeNavigation(html, "Help", "help/index.html");
    }

    html << "<div style='clear:both;height:10px;'></div>";

}


void ExpressionMatrix::writeNavigation(ostream& html, const string& text, const string& url)
{

    html <<
        "<button onclick=\"window.location = '" << url << "';\""
        " style='"
        "background-color:LightCoral;"
        "width:100px;"
        "height:50px;"
        "margin:3px;"
        "border-radius:6px;"
        "text-align:center;"
        "vertical-align:middle;"
        "float:left;"
        "cursor:pointer;"
        "'>\n"
        << text << "</button>";
}



void ExpressionMatrix::exploreSummary(
    const vector<string>& request,
    ostream& html)
{



    // If the run directory contains a README file, copy it to html.
    // The README file can contain html, if it is well behaved.
    if(boost::filesystem::exists("README") && boost::filesystem::is_regular_file("README")) {
    	ifstream readMeFile("README");
    	html << "<h1>README file</h1>";
    	html << readMeFile.rdbuf();
    }



    html <<
        "<h1>Expression matrix</h1>"
        "<p>The expression matrix has " << cellCount() <<
        " cells and " << geneCount() <<
        " genes.";


    // Link to the summary of hash table usage.
    html << "<p>See a <a href=exploreHashTableSummary>summary of hash table utilization</a>.";
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
        const bool isBad = (loadFactor >= 0.5);
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
        // Hash table used to hold cell names.
        const double capacity = double(cellNames.capacity());
        const double size = double(cellCount());
        const double loadFactor = size / capacity;
        const double penalty = 0.5 * (1. + 1. / ((1.-loadFactor) * (1.-loadFactor))) - 1.;
        const bool isBad = (loadFactor >= 0.5);
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
        const bool isBad = (loadFactor >= 0.5);
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
        const bool isBad = (loadFactor >= 0.5);
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



void ExpressionMatrix::exploreCell(
    const vector<string>& request,
    ostream& html)
{

    // Get the cell id.
    string cellIdString;
    const bool cellIdIsPresent = getParameterValue(request, "cellId", cellIdString);
    const CellId cellId = cellIdFromString(cellIdString);

    // Write the form to get the cell id.
    html <<
        "<form>"
        "Specify a cell using a case-sensitive name or a numeric cell id between 0 and " << cellCount()-1 <<
        " included:<br><input type=text name=cellId autofocus>"
        "</form>";

    // If there is no cell id, do nothing.
    if(!cellIdIsPresent) {
        return;
    }

    // Access the cell.
    if(cellId==invalidCellId) {
        html << "<p>Invalid cell id";
        return;
    }
    const Cell& cell = cells[cellId];
    const string& cellName = cellNames[cellId];

    // Write a title.
    html << "<h1>Cell " << cellId << " " << cellName << "</h1>";


    // Write a table containing meta data and additional information for this cell.
    html << "<h2>Cell meta data and additional cell information</h2>";
    html << "<p><table>";
    for(const auto& p: cellMetaData[cellId]) {
        html << "<tr><td>" << cellMetaDataNames[p.first] << "<td>" << cellMetaDataValues[p.second];
    }
    html << "<tr><td>Cell id<td>" << cellId;
    html << "<tr><td>Number of genes with non-zero expression counts<td>" <<
        cellExpressionCounts.size(cellId);
    html << "<tr><td>Sum of expression counts<td>" << cell.sum1;
    html << "</table>";

    // The expression counts are stored sorted by gene id,
    // but we want to show them in order of decreasing count.
    const auto storedExpressionCounts = cellExpressionCounts[cellId];
    vector< pair<GeneId, float> > expressionCounts(
        storedExpressionCounts.begin(),
        storedExpressionCounts.end());
    sort(
        expressionCounts.begin(),
        expressionCounts.end(),
        OrderPairsBySecondGreaterThenByFirstLess< pair<GeneId, float> >());



    // Write to html jQuery and TableSorter so we can make the table below sortable.
    writeJQuery( html);
    writeTableSorter(html);



    // Write a table of the expression counts for this cell.
    html << "<h2>Gene expression counts for this cell</h2>";
    html <<
        "<p>The following table of expression counts for this cell is sortable. Click on a header to sort by that header. "
        "Click again to reverse the sorting order."
        "<p><table id=countTable class=tablesorter><thead><tr><th>Gene<br>name<th>Raw<br>count"
            "<th>L1-normalized<br>count<br>(sum is 1)"
            "<th>L2-normalized<br>count<br>(sum<br>of<br>squares is 1)</thead><tbody>";
    for(const auto& p: expressionCounts) {
        const GeneId geneId = p.first;
        CZI_ASSERT(geneId < geneCount());
        const string geneName = geneNames[geneId];
        const float count = p.second;
        html <<  "<tr><td class=centered><a href=gene?geneId=" << urlEncode(geneName) << ">" << geneName << "</a>";
        html <<
            "<td class=centered>" << count;
        const auto oldPrecision = html.precision(3);
        html <<
            "<td class=centered>" << count * cell.norm1Inverse <<
            "<td class=centered>" << count * cell.norm2Inverse;
        html.precision(oldPrecision);
    }

    // Finish the table and make it sortable.
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#countTable').tablesorter();});"
        "</script>"
        ;

}



void ExpressionMatrix::compareTwoCells(
    const vector<string>& request,
    ostream& html)
{
    // Get the cell ids.
    array<string, 2> cellIdStrings;
    const bool cellId0IsPresent = getParameterValue(request, "cellId0", cellIdStrings[0]);
    const bool cellId1IsPresent = getParameterValue(request, "cellId1", cellIdStrings[1]);
    const bool cellIdsArePresent = cellId0IsPresent && cellId1IsPresent;

    // Write the form to get the cell ids.
    html <<
        "<form>"
        "Specify two cells using names or numeric ids between 0 and " << cells.size()-1 << ":"
        "<br><input type=text name=cellId0 autofocus";
    if(cellId0IsPresent) {
        html << " value=" << cellIdStrings[0];
    }
    html <<
        ">"
        "<br><input type=text name=cellId1";
    if(cellId1IsPresent) {
        html << " value=" << cellIdStrings[1];
    }
    html <<
        ">"
        "<input type=submit hidden>"
        "</form>";

    // If the cell ids are not specified, do nothing.
    if(!cellIdsArePresent) {
        return;
    }

    // Access the cells.
    array<CellId, 2> cellIds;
    for(int i=0; i<2; i++) {
        cellIds[i] = cellIdFromString(cellIdStrings[i]);
        if(cellIds[i]<0) {
            html << "<p>Invalid cell id " << cellIdStrings[i];
            return;
        }
    }
    // const Cell& cell0 = cells[cellIds[0]];
    // const Cell& cell1 = cells[cellIds[1]];

    // Write a title.
    html  << "<h1>Comparison of cells " << cellIds[0] << " ";
    writeCellLink(html, cellIds[0], false);
    html << " and " << cellIds[1] << " ";
    writeCellLink(html, cellIds[1], false);
    html << "</h1>";


    // Write a table of similarities between these two cells.
    html << "<table>";
    html << "<tr><th class=left>Exact similarity<td>" << computeCellSimilarity(cellIds[0], cellIds[1]);
    html << "<tr><th class=left>Approximate similarity<td>" << computeApproximateCellSimilarity(cellIds[0], cellIds[1]);
    html << "</table>";



    // Create a table of (totalCount, geneId, count for cell0, count for cell1).
    vector< tuple<double, int, double, double> > data;
    const auto& cellCounts0 = cellExpressionCounts[cellIds[0]];
    const auto& cellCounts1 = cellExpressionCounts[cellIds[1]];
    const auto begin0 = cellCounts0.begin();
    const auto begin1 = cellCounts1.begin();
    const auto end0   = cellCounts0.end();
    const auto end1   = cellCounts1.end();
    auto it0 = begin0;
    auto it1 = begin1;
    while(true) {
        if(it0==end0) {
            if(it1==end1) {
                // Both it0 and it1 are at end. We are done.
                break;
            } else {
                // it0 is at end, it1 is not.
                const size_t geneId = it1->first;
                const double count1 = it1->second;
                data.push_back(make_tuple(count1, geneId, 0, count1));
                ++it1;
            }
        } else {
            if(it1==end1) {
                // it1 is at end, it0 is not.
                const size_t geneId = it0->first;
                const double count0 = it0->second;
                data.push_back(make_tuple(count0, geneId, count0, 0));
                ++it0;
            } else {
                // Neither it0 nor it1 are at end.
                if(it0->first < it1->first) {
                    const size_t geneId = it0->first;
                    const double count0 = it0->second;
                    data.push_back(make_tuple(count0, geneId, count0, 0));
                    ++it0;
                } else if(it1->first < it0->first) {
                    const size_t geneId = it1->first;
                    const double count1 = it1->second;
                    data.push_back(make_tuple(count1, geneId, 0, count1));
                    ++it1;
                } else {
                    const size_t geneId = it0->first;
                    const double count0 = it0->second;
                    const double count1 = it1->second;
                    data.push_back(make_tuple(count0+count1, geneId, count0, count1));
                    ++it0;
                    ++it1;
                }
            }
        }
    }

    // Sort by decreasing total count.
    sort(data.begin(), data.end(), std::greater< tuple<double, int, double, double> >());


    // Compute the maximum counts.
    double maxCount0 = 0;
    double maxCount1 = 0;
    for(const auto& t: data) {
        const double count0 = t.get<2>();
        const double count1 = t.get<3>();
        maxCount0 = max(maxCount0, count0);
        maxCount1 = max(maxCount1, count1);
    }





    // Draw a scatter plot of the expression counts for the two cells.
    html <<
        "<script src='https://www.gstatic.com/charts/loader.js'></script>"
        "<script>"
        "    google.charts.load('current', {'packages':['corechart']});"
        "    google.charts.setOnLoadCallback(drawChart);"

        "    var xMax = " << maxCount0 << ";"
        "    var yMax = " << maxCount1 << ";"

        "    function roundUp(x)"
        "    {"
        "        var y = Math.pow(10., Math.floor(Math.log10(x)));"
        "        if(x <=2*y) {"
        "            return 2*y;"
        "        } else if(x<=5*y) {"
        "            return 5*y;"
        "        } else {"
        "            return 10*y;"
        "        }"
        "    }"

        "    var data;"
        "    var chart;"
        "    var options = {"
        "        hAxis: {title: 'Count for cell " << cellIds[0] << " " << cellNames[cellIds[0]] << "', viewWindowMode: 'explicit', viewWindow: {min:0, max:0}},"
        "        vAxis: {title: 'Count for cell " << cellIds[1] << " " << cellNames[cellIds[1]] << "', viewWindowMode: 'explicit', viewWindow: {min:0, max:0}},"
        "        legend: 'none',"
        "        pointSize: 2,"
        "        tooltip: {isHtml: true, trigger: 'both'}"
        "    };"


        "    function drawChart()"
        "    {"
        "        data = new google.visualization.DataTable();"
        "        data.addColumn('number', 'x');"
        "        data.addColumn('number', 'y');"
        "        data.addColumn({type: 'string', role: 'tooltip', 'p': {'html': true}});"
        "        data.addRows([";
    for(const auto& t: data) {
        const int geneId = t.get<1>();
        const double count0 = t.get<2>();
        const double count1 = t.get<3>();
        html << "[" << count0 << "," << count1 << ",";
        // html << "<span onclick=\"window.location = \\x22gene?geneId=" << geneId << "\\x22\">" << genes[geneId].name << "</span>";
        html << "\"<a href='gene?geneId=" << geneId << "'>" << geneNames[geneId] << "</a>\"";
        html << "],";
    }
    html <<
        "        ]);"
        "        chart = new google.visualization.ScatterChart(document.getElementById('scatterPlot'));"
        "        options.hAxis.viewWindow.max = roundUp(xMax);"
        "        options.vAxis.viewWindow.max = roundUp(yMax);"
        "        chart.draw(data, options);"
        "    }"

        "    function scale(factor)"
        "    {"
        "        xMax = factor * xMax;"
        "        yMax = factor * yMax;"
        "        options.hAxis.viewWindow.max = roundUp(xMax);"
        "        options.vAxis.viewWindow.max = roundUp(yMax);"
        "        chart.draw(data, options);"
        "    }"

        "    function horizontalScale(factor)"
        "    {"
        "        xMax = factor * xMax;"
        "        options.hAxis.viewWindow.max = roundUp(xMax);"
        "        chart.draw(data, options);"
        "    }"

        "    function verticalScale(factor)"
        "    {"
        "        yMax = factor * yMax;"
        "        options.vAxis.viewWindow.max = roundUp(yMax);"
        "        chart.draw(data, options);"
        "    }"

        "    function changeMarkerSize(increment)"
        "    {"
        "        options.pointSize += increment;"
        "        if(options.pointSize<1) {"
        "            options.pointSize = 1.;"
        "        }"
        "        chart.draw(data, options);"
        "    }"

        "</script>"
        "<div id='scatterPlot' style='float: left;width: 800px; height: 600px;'></div>"
        "<div style='float:left'>"
        "<input type=button value='Zoom in' onclick='scale(0.5);' style='width:180px;border-radius:5px;' /><br>"
        "<input type=button value='Zoom out' onclick='scale(2);' style='width:180px;border-radius:5px;' /><br>"
        "<input type=button value='Horizontal only zoom in' onclick='horizontalScale(0.5);' style='width:180px;border-radius:5px;' /><br>"
        "<input type=button value='Horizontal only zoom out' onclick='horizontalScale(2);' style='width:180px;border-radius:5px;' /><br>"
        "<input type=button value='Vertical only zoom in' onclick='verticalScale(0.5);' style='width:180px;border-radius:5px;' /><br>"
        "<input type=button value='Vertical only zoom out' onclick='verticalScale(2);' style='width:180px;border-radius:5px;' /><br>"
        "<input type=button value='Larger marker' onclick='changeMarkerSize(1);' style='width:180px;border-radius:5px;' /><br>"
        "<input type=button value='Smaller marker' onclick='changeMarkerSize(-1);' style='width:180px;border-radius:5px;' /><br>"
        "</div>"
        "<div style='clear:both;' />"
        ;



    // Write the table with the counts.
    html <<
        "<table>"
        "<tr><th>Gene<br>id<th>Gene<br>name<th>Count for<br>cell " << cellIds[0] << "<br>";
    writeCellLink(html, cellIds[0], false);
    html << "<th>Count for<br>cell " << cellIds[1] << "<br>";
    writeCellLink(html, cellIds[1], false);

    for(const auto& t: data) {
        const int geneId = t.get<1>();
        const double count0 = t.get<2>();
        const double count1 = t.get<3>();
        html << "<tr><td>";
        writeGeneLink(html, geneId, true);
        html << "<td>";
        writeGeneLink(html, geneId, false);
        html << "<td>" << count0 << "<td>" << count1;
    }
    html << "</table>";


}



// Write a hyperlink for a cell. The last parameter controls whether the link is
// written as a cell id or cell name.
ostream& ExpressionMatrix::writeCellLink(ostream& html, CellId cellId, bool writeId)
{
    if(cellId == invalidCellId) {
        html << "Invalid cell";
    } else {
        CZI_ASSERT(cellId < cells.size());
        html << "<a href='cell?cellId=" << cellId << "'>";
        if(writeId) {
            html << cellId;
        } else {
            html << cellNames[cellId];
        }
        html << "</a>";
    }
    return html;
}
ostream& ExpressionMatrix::writeCellLink(ostream& html, const string& cellName, bool writeId)
{
    const CellId cellId = cellIdFromString(cellName);
    writeCellLink(html, cellId, writeId);
    return html;
}



// Write a hyperlink for a gene. The last parameter controls whether the link is
// written as a gene id or gene name.
ostream& ExpressionMatrix::writeGeneLink(ostream& html, GeneId geneId, bool writeId)
{
    if(geneId == invalidGeneId) {
        html << "Invalid gene";
    } else {
        CZI_ASSERT(geneId < geneNames.size());
        html << "<a href='gene?geneId=" << geneId << "'>";
        if(writeId) {
            html << geneId;
        } else {
            html << geneNames[geneId];
        }
        html << "</a>";
    }
    return html;
}
ostream& ExpressionMatrix::writeGeneLink(ostream& html, const string& geneName, bool writeId)
{
    const GeneId geneId = geneIdFromString(geneName);
    writeGeneLink(html, geneId, writeId);
    return html;
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
    html << " style='width:8em;'>";


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



ostream& ExpressionMatrix::writeCellSetSelection(
    ostream& html,
    const string& selectName,
    bool multiple) const
{
    set<string> selected;
    return writeCellSetSelection(html, selectName, selected, multiple);
}

ostream& ExpressionMatrix::writeCellSetSelection(
    ostream& html,
    const string& selectName,
    const set<string>& selected,
    bool multiple) const
{
    html << "<select";
    if(multiple) {
        html << " multiple";
    }
    html << " name=" << selectName << " style='vertical-align:text-top;'>";
    html << "<option value=''></option>";
    for(const auto& p: cellSets.cellSets) {
        const string& cellSetName = p.first;
        html << "<option value=" << cellSetName;
        if(selected.find(cellSetName) != selected.end()) {
            html << " selected=selected";
        }
        html << ">" << cellSetName << "</option>";
    }
    html << "</select>";

    return html;
}



ostream& ExpressionMatrix::writeGraphSelection(
    ostream& html,
    const string& selectName,
    bool multiple) const
{
    html << "<select";
    if(multiple) {
        html << " multiple";
    }
    html << " name=" << selectName << "><option value=''></option>";
    for(const auto& p: graphs) {
        const string& graphName = p.first;
        html << "<option value=" << graphName << ">" << graphName << "</option>";
    }
    html << "</select>";

    return html;
}



void ExpressionMatrix::exploreCellSets(
    const vector<string>& request,
    ostream& html)
{
    // Write a title.
    html << "<h1>Cell sets</h1>";

    // Write a table listing the cell sets in existence.
    html << "<p><table><th>Cell<br>set<br>name<th>Number<br>of<br>cells<th class=centered>Click<br>to<br>remove";
    for(const auto& p: cellSets.cellSets) {
        const string& name = p.first;
        const auto& cellSet = *p.second;
        html << "<tr><td><a href='cellSet?cellSetName=" << urlEncode(name) << "'>" << name << "</a><td class=centered>" << cellSet.size();
        html << "<td  class=centered>";
        if(name != "AllCells") {
            html << "<a href='removeCellSet?cellSetName=" << urlEncode(name) << "'>Remove</a>";
        }
    }
    html << "</table>";


    // Form to create a new cell set from meta data.
    html <<
        "<br><h2>Create a new cell set using meta data</h2>"
        "<p><form action=createCellSetUsingMetaData>"
        "<input type=submit value='Create a new cell set'> with name "
        "<input type=text required name=cellSetName>"
        " consisting of cells for which meta data field ";
        set<string> metaDataNames;
        writeMetaDataSelection(html, "metaData", metaDataNames, false);
    html <<
        " matches this regular expression: "
        "<input type=text name=regex>"
        "</form>";




    // Form to create a new cell set by union/intersection existing cell sets.
    html <<
        "<br><h2>Create a new cell set by union/intersection of existing cell sets</h2>"
        "<p><form action=createCellSetIntersectionOrUnion>"
        "<input type=submit value='Create a new cell set'> with name "
        "<input type=text required name=cellSetName>"
        " as the "
        "<select name=operation>"
        "<option value=union>union</option>"
        "<option value=intersection>intersection</option>"
        "</select>"
        " of the selected cell sets: ";
    writeCellSetSelection(html, "inputCellSets", true);
    html << "</form>";



    // Form to create a new cell set by downsampling an existing cell set.
    html <<
        "<br><h2>Create a new cell set by downsampling an existing cell set</h2>"
        "<p>The new cell set will be a random subset of the specified cell set."
    	" Each cell in the specified cell set is inserted in the random subset with the specified probability."
    	" Therefore, the downsampling rate will be approximately equal to the specified probability."
        "<p><form action=downsampleCellSet>"
        "<input type=submit value='Create a new cell set'> with name "
        "<input type=text required name=cellSetName>"
        " by downsampling cell set ";
    writeCellSetSelection(html, "inputCellSet", false);
    html <<
		" with probability "
		"<input type=text required name=probability size=6>"
		" and random seed "
		"<input type=text required name=seed value='231' size=6>"
    	"</form>";

}



void ExpressionMatrix::exploreCellSet(
    const vector<string>& request,
    ostream& html)
{
    // Get the name of the cell set we want to look at.
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "Missing cell set name.";
        return;
    }

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

    // Write a title.
    html << "<h1>Cell set " << cellSetName << "</h1>";

    // Locate the cell set.
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        html << "<p>This cell set does not exist.";
        return;
    }
    const auto& cellSet = *it->second;
    html << "<p>This cell set has " << cellSet.size() << " cells." << endl;



    // Write the form to get the metadata to display.
    html <<
        "<form>"
        "Select cell metadata to display:<br>";
    writeMetaDataSelection(html, "metadata", metaDataToDisplay, true);
    html <<
        "<input type=hidden name=cellSetName value='" << cellSetName << "'>"
        "<br><input type=submit value='Redisplay table'>"
        "</form>";



    // Write a table containing the cells of this set.
    html << "<br><table><tr><th class=centered>Cell<br>id<th class=centered>Cell<br>name";
    for(const auto& metaDataFieldName: metaDataToDisplayStrings) {
        html << "<th>" << metaDataFieldName.second;
    }
    for(const CellId cellId: cellSet) {
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
    html << "</table>";
}



void ExpressionMatrix::createCellSetUsingMetaData(const vector<string>& request, ostream& html)
{
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "Missing cell set name.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }

    string metaData;
    if(!getParameterValue(request, "metaData", metaData)) {
        html << "Missing meta data name.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }


    string regex;
    if(!getParameterValue(request, "regex", regex)) {
        html << "Missing regular expression.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }
    string decodedRegex;
    urlDecode(regex, decodedRegex);

    if(createCellSetUsingMetaData(cellSetName, metaData, decodedRegex)) {
        html << "<p>Newly created cell set " << cellSetName << " has ";
        html << cellSets.cellSets[cellSetName]->size() << " cells.";
    } else {
        html << "<p>Unable to create cell set " << cellSetName << ".";
    }
    html << "<p><form action=cellSets><input type=submit value=Continue></form>";

}



void ExpressionMatrix::createCellSetIntersectionOrUnion(const vector<string>& request, ostream& html)
{
    // Get the name of the cell set to be created.
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "Missing cell set name.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }

    // Get the name of the operation to be performed (intersection or union).
    string operation;
    if(!getParameterValue(request, "operation", operation)) {
        html << "Missing operation.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }
    bool doUnion;
    if(operation == "intersection") {
        doUnion = false;
    } else if(operation == "union") {
        doUnion = true;
    } else {
        html << "Invalid operation.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }


    // Get the names of the input cell sets.
    set<string> inputCellSets;
    getParameterValues(request, string("inputCellSets"), inputCellSets);
    if(inputCellSets.size() < 2) {
        html << "At least two input cell sets should be specified.";
        return;
    }

    // Concatenate the input cell sets with commas.
    string inputCellSetsString;
    for(const string& inputCellSet: inputCellSets) {
        inputCellSetsString.append(inputCellSet);
        inputCellSetsString.append(",");
    }
    inputCellSetsString.resize(inputCellSetsString.size()-1);


    // Do the intersection or union.
    if(createCellSetIntersectionOrUnion(inputCellSetsString, cellSetName, doUnion)) {
        html << "<p>Newly created cell set " << cellSetName << " has ";
        html << cellSets.cellSets[cellSetName]->size() << " cells.";
    } else {
        html << "<p>Unable to create cell set " << cellSetName << ".";
    }
    html << "<p><form action=cellSets><input type=submit value=Continue></form>";

}



void ExpressionMatrix::downsampleCellSet(const vector<string>& request, ostream& html)
{
    // Get the name of the cell set to be created.
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "Missing cell set name.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }


    // Get the name of the input cell set.
    string inputCellSet;
    getParameterValue(request, "inputCellSet", inputCellSet);

    // Get the downsampling parameters.
    double probability;
    getParameterValue(request, "probability", probability);
    int seed;
    getParameterValue(request, "seed", seed);


    // Do the downsampling.
    if(downsampleCellSet(inputCellSet, cellSetName, probability, seed)) {
        html << "<p>Newly created cell set " << cellSetName << " has ";
        html << cellSets.cellSets[cellSetName]->size() << " cells.";
        html << "<p>Downsampling probability was " << probability;
        html << "<p>Actual downsampling rate was " << double(cellSets.cellSets[cellSetName]->size()) / double(cellSets.cellSets[inputCellSet]->size());
    } else {
        html << "<p>Unable to create cell set " << cellSetName << ".";
    }


    // The button to continue goes back to the cell sets page.
    html << "<p><form action=cellSets><input type=submit value=Continue></form>";

}



void ExpressionMatrix::removeCellSet(const vector<string>& request, ostream& html)
{
    // Get the name of the cell set we want to look at.
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "Missing cell set name.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }

    // Locate the cell set.
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        html << "<p>Cell set " << cellSetName << " cannot be removed because it does not exist.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }


    const string fileName = it->second->fileName;
    cellSets.cellSets.erase(it);
    if(boost::filesystem::remove(fileName)) {
        html << "<p>Cell set " << cellSetName << " was removed.";
    } else {
        html << "<p>Cell set " << cellSetName << " was removed from memory but the corresponding memory mapped file ";
        html << fileName << " could not be removed.";
    }
    html << "<p><form action=cellSets><input type=submit value=Continue></form>";

}



void ExpressionMatrix::exploreGene(
    const vector<string>& request,
    ostream& html)
{

    // Get the gene id.
    string geneIdString;
    const bool geneIdIsPresent = getParameterValue(request, "geneId", geneIdString);
    const GeneId geneId = geneIdFromString(geneIdString);

    // Get the name of the cell set to use.
    // We will only write expression counts for cells in this cell set.
    string cellSetName;
    getParameterValue(request, "cellSetName", cellSetName);

    // Get the names and string ids of the meta data to display.
    vector<string> metaDataToDisplay;
    getParameterValues(request, string("metaDataName"), metaDataToDisplay);
    vector<StringId> metaDataToDisplayStringIds;
    for(const string& metaDataName: metaDataToDisplay) {
        const StringId stringId = cellMetaDataNames(metaDataName);
        if(stringId != cellMetaDataNames.invalidStringId) {
            metaDataToDisplayStringIds.push_back(cellMetaDataNames(metaDataName));
        }
    }

    // Write the form.
    html <<
        "<form>"
        "Specify a gene using a case-sensitive name or a numeric id between 0 and  " << geneCount()-1 <<
        " included:<br><input type=text name=geneId";
    if(geneIdIsPresent) {
        html << " value=" << geneIdString;
    }
    html << " autofocus>";
    html << "<br>Write expression counts for cells in cell set ";
    const set<string> selectedCellSet = {cellSetName};
    writeCellSetSelection(html, "cellSetName", selectedCellSet, false);
    html << "<br>Display these cell meta data fields:";
    writeMetaDataSelection(html, "metaDataName", metaDataToDisplay, true);
    html << "<br><input type=submit value='Submit'>";
    html << "</form>";



    // If there is no gene id, do nothing.
    if(!geneIdIsPresent) {
        return;
    }

    // Access the gene.
    if(geneId == invalidGeneId) {
        html << "<p>Invalid gene id";
        return;
    }
    const string& geneName = geneNames[geneId];


    // Write a title.
    html << "<h1>Gene "<< geneId << " " << geneName << "</h1>";

    // Links with information on this gene on various web sites.
    html <<
        "<p>Look up this gene on "
        "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=" << urlEncode(geneName) << "'>NCBI</a>, "
        "<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg38&position=" << urlEncode(geneName) << "'>UCSC</a>, "
        "<a href='http://www.ensembl.org/Multi/Search/Results?q=" << urlEncode(geneName) << ";site=ensembl'>Ensembl</a>, "
        "<a href='https://en.wikipedia.org/wiki/" << urlEncode(geneName) << "'>Wikipedia</a>, "
        "<a href='https://www.google.com/#q=gene+" << urlEncode(geneName) << "'>Google</a>."
        ;

    // If there is no cell set name, stop here.
    if(cellSetName.empty()) {
        html << "<p>To get expression counts for this gene, specify a cell set in the form above.";
        return;
    }

    // Locate the cell set.
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        html << "Cell set " << cellSetName << " does not exists.";
        return;
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);






    // Gather the expression counts for this gene for all cells in the specified cell set.
    vector<ExploreGeneData> counts;
    for(const CellId cellId: cellSet) {
        const Cell& cell = cells[cellId];
        ExploreGeneData data;
        data.cellId = cellId;
        data.rawCount =  getExpressionCount(cellId, geneId);
        if(data.rawCount) {
            data.count1 = float(data.rawCount * cell.norm1Inverse);
            data.count2 = float(data.rawCount * cell.norm2Inverse);
            counts.push_back(data);
        }
    }


    // Sort the counts.
    // For now this just sorts by L2-normalized count. but later
    // we should add options for sorting in other ways.
    sort(counts.begin(), counts.end());

    // Write to html jQuery and TableSorter so we can make the table below sortable.
    writeJQuery( html);
    writeTableSorter(html);



    // Write a table with the counts.
    // WE ALSO HAVE TO ADD THE OUTPUT OF THE REQUESTED META DATA.
    html <<
        "<p>The following table of expression counts for this gene is sortable. Click on a header to sort by that header. "
        "Click again to reverse the sorting order."
        "<br><table id=countTable class=tablesorter><thead><tr><th>Cell<br>id<th>Cell<br>name"
        "<th>Unnormalized<br>count<br>(raw<br>count)"
        "<th>L1-normalized<br>count<br>(fractional<br>read<br>count)"
        "<th>L2-normalized<br>count";
    for(const StringId metaDataNameStringId: metaDataToDisplayStringIds) {
        html << "<th>" << cellMetaDataNames[metaDataNameStringId];
    }
    html << "</thead><tbody>";
    for(const auto& data: counts) {
        html << "<tr><td class=centered>";
        writeCellLink(html, data.cellId, true);
        html <<"<td>";
        writeCellLink(html, data.cellId, false);
        html << "<td class=centered>" << data.rawCount;
        const auto oldPrecision = html.precision(3);
        html << "<td class=centered>" << data.count1;
        html << "<td class=centered>" << data.count2;
        html.precision(oldPrecision);

        // Write the requested meta data for this cell.
        for(const StringId metaDataNameStringId: metaDataToDisplayStringIds) {
            html << "<td class=centered>" << getMetaData(data.cellId, metaDataNameStringId);
        }
    }


    // Finish the table and make it sortable.
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#countTable').tablesorter();});"
        "</script>"
        ;
}



void ExpressionMatrix::exploreGraphs(
    const vector<string>& request,
    ostream& html)
{
    html << "<h1>Cell similarity graphs</h1>";



    // Table with existing graphs.
    html <<
        "<table><tr>"
        "<th class=centered>Graph<br>name"
        "<th class=centered>Cell<br>set<br>name"
        "<th class=centered>Similar<br>pairs<br>name"
        "<th class=centered>Similarity<br>threshold"
        "<th class=centered>Maximum<br>connectivity"
        "<th class=centered>Number<br>of<br>vertices"
		"<th class=centered>Number<br>of<br>edges"
		"<th class=centered>Number<br>of<br>isolated<br>vertices<br>removed"
        "<th class=centered>Clustering"
        "<th class=centered>Action";
    for(const auto& p: graphs) {
        const string& graphName = p.first;
        const GraphInformation& info = p.second.first;
        // const CellSimilarityGraph& graph = *(p.second.second);
        html << "<tr><td><a href='graph?graphName=" << urlEncode(graphName);
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
        html << "<td class=centered><form action=clusterDialog><input type=text hidden name=graphName value='" << graphName << "'><input type=submit value='Run clustering on graph " << graphName << "'></form>";
        html << "<td class=centered><form action=removeGraph><input type=text hidden name=graphName value='" << graphName << "'><input type=submit value='Remove graph " << graphName << "'></form>";
    }



    // Add a final row to the table, containing a form to create a new graph.
    vector<string> availableSimilarPairs;
    getAvailableSimilarPairs(availableSimilarPairs);

    html <<
        "<tr>"
        "<form action=createNewGraph>"
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
        "<td><td><td><td><td class=centered><input type=submit value='Create a new graph'>"
        "</form>";

    html << "</table>";



    // Form to compare two graphs.
    if(graphs.size() > 1) {
        html << "<p><form action=compareGraphs><input type=submit value=Compare> graphs ";
        writeGraphSelection(html, "graphName0", false);
        html << " and ";
        writeGraphSelection(html, "graphName1", false);
        html << ".</form>";
    }




}



void ExpressionMatrix::compareGraphs(
    const vector<string>& request,
    ostream& html)
{
    // Get the names of the two graphs to be compared.
    string graphName0;
    getParameterValue(request, "graphName0", graphName0);
    string graphName1;
    getParameterValue(request, "graphName1", graphName1);

    // Locate the graphs.
    const auto it0 = graphs.find(graphName0);
    const auto it1 = graphs.find(graphName1);
    if(it0==graphs.end() || it1==graphs.end()) {
        html << "<p>Did not find one or both of the graphs to be compared.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }
    if(it0 == it1) {
        html << "<p>Comparison of a graph with itself requested. ";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }
    const GraphInformation& graphCreationParameters0 = it0->second.first;
    const GraphInformation& graphCreationParameters1 = it1->second.first;
    const auto graphPointer0 = it0->second.second;
    const auto graphPointer1 = it1->second.second;
    CZI_ASSERT(graphPointer0);
    CZI_ASSERT(graphPointer1);
    const CellSimilarityGraph& graph0 = *graphPointer0;
    const CellSimilarityGraph& graph1 = *graphPointer1;



    // Find the common vertices (vertices that correspond to the same cell).
    vector<CellId> cells0, cells1;
    BGL_FORALL_VERTICES(v, graph0, CellSimilarityGraph) {
        cells0.push_back(graph0[v].cellId);
    }
    BGL_FORALL_VERTICES(v, graph1, CellSimilarityGraph) {
        cells1.push_back(graph1[v].cellId);
    }
    sort(cells0.begin(), cells0.end());
    sort(cells1.begin(), cells1.end());
    vector<CellId> commonCells;
    std::set_intersection(cells0.begin(), cells0.end(), cells1.begin(), cells1.end(), back_inserter(commonCells));



    // Maps of the edges. Keyed by pair(CellId, CellId), with the lowest numbered cell first.
    // Values: similarities.
    map< pair<CellId, CellId>, float> edgeMap0, edgeMap1;
    BGL_FORALL_EDGES(e, graph0, CellSimilarityGraph) {
        const CellSimilarityGraph::vertex_descriptor vA = source(e, graph0);
        const CellSimilarityGraph::vertex_descriptor vB = target(e, graph0);
        CellId cellIdA = graph0[vA].cellId;
        CellId cellIdB = graph0[vB].cellId;
        CZI_ASSERT(cellIdA != cellIdB);
        if(cellIdB < cellIdA) {
            swap(cellIdA, cellIdB);
        }
        CZI_ASSERT(cellIdA < cellIdB);
        edgeMap0.insert(make_pair( make_pair(cellIdA, cellIdB), graph0[e].similarity));
    }
    BGL_FORALL_EDGES(e, graph1, CellSimilarityGraph) {
        const CellSimilarityGraph::vertex_descriptor vA = source(e, graph1);
        const CellSimilarityGraph::vertex_descriptor vB = target(e, graph1);
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




void ExpressionMatrix::exploreGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the graph name.
    string graphName;
    if(!getParameterValue(request, "graphName", graphName)) {
        html << "Missing graph name.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }

    // Find the graph.
    const auto it = graphs.find(graphName);
    if(it == graphs.end()) {
        html << "<p>Graph " << graphName << " does not exists.";
        return;
    }
    const GraphInformation& graphInfo = it->second.first;
    CellSimilarityGraph& graph = *(it->second.second);

    // Write the title.
    html << "<h1>Graph " << graphName << "</h1>";


    // Compute the graph layout, if necessary.
    if(!graph.layoutWasComputed) {
        html << "<div style='font-family:courier'>";
		html << timestamp << "Graph layout computation begins.";
		graph.computeLayout();
		html << "<br>" << timestamp << "Graph layout computation ends.";
		html << "</div>";
		graph.layoutWasComputed = true;
    }


    // Div to contain the table and the form for coloring.
    html << "<div>";

    // Write a table with the graph creation parameters.
    html << "<div style='float:left;margin:10px'>";
    html << "<table>";
    html << "<tr><td>Cell set name<td><a href='cellSet?cellSetName=" << urlEncode(graphInfo.cellSetName);
    html << "'>" << graphInfo.cellSetName << "</a>";
    html << "<tr><td>Similar pairs name<td>" << graphInfo.similarPairsName;
    html << "<tr><td>Similarity threshold<td class=centered>" << graphInfo.similarityThreshold;
    html << "<tr><td>Maximum connectivity<td class=centered>" << graphInfo.maxConnectivity;
    html << "<tr><td>Number of vertices (cells)<td class=centered>" << boost::num_vertices(graph);
    html << "<tr><td>Number of edges<td class=centered>" << boost::num_edges(graph);
    html << "<tr><td>Number of isolated vertices (cells) removed<td class=centered>" << graphInfo.isolatedVertexCount;
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
    string normalizationOption;
    getParameterValue(request, "normalizationOption", normalizationOption);
    if(normalizationOption.empty()) {
        normalizationOption = "norm2";
    }
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
        "    normalizationOptionChanged = document.getElementById('normalizationOption').value != '" << normalizationOption << "';\n"
        "    // window.alert('Normalization option changed: ' + normalizationOptionChanged);\n"
        "    cellIdChanged = document.getElementById('cellIdInput').value != '" << cellIdStringForColoringBySimilarity << "';\n"
        "    // window.alert('Cell id changed: ' + cellIdChanged);\n"
        "    metaDataChanged = document.getElementById('metaDataName').value != '" << metaDataName << "';\n"
        "    // window.alert('Meta data changed: ' + metaDataChanged);\n"
        "    if(coloringOptionChanged || geneIdChanged || normalizationOptionChanged || cellIdChanged || metaDataChanged) {\n"
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

    html << " using <select name=normalizationOption id=normalizationOption>";
    html << "<option value=none";
    if(normalizationOption == "none") {
        html << " selected=selected";
    }
    html << ">no normalization (raw read count)</option>";
    html << "<option value=norm1";
    if(normalizationOption == "norm1") {
        html << " selected=selected";
    }
    html << ">L1 normalization (fractional read count)</option>";
    html << "<option value=norm2";
    if(normalizationOption == "norm2") {
        html << " selected=selected";
    }
    html << ">L2 normalization</option>";
    html << "</select>";




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
&#129094;
</button>
<button onClick='moveLeft()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
&#129092;
</button>
<button onClick='moveUp()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
&#129093;
</button>
<button onClick='moveDown()' style='width:120px;height:30px;margin:2px;vertical-align:middle;border-radius:6px;background-color:pink;'>
&#129095;
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
        BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
            CellSimilarityGraphVertex& vertex = graph[v];
            const double rawCount = getExpressionCount(vertex.cellId, geneId);
            if(normalizationOption == "none") {
                vertex.value = rawCount;
            } else {
                const Cell& cell = cells[vertex.cellId];
                if(normalizationOption == "norm1") {
                    vertex.value = rawCount * cell.norm1Inverse;
                } else if(normalizationOption == "norm2") {
                    vertex.value = rawCount * cell.norm2Inverse;
                } else {
                    html << "<p>Invalid normalization option.";
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
        BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
            CellSimilarityGraphVertex& vertex = graph[v];
            vertex.value = computeCellSimilarity(cellIdForColoringBySimilarity, vertex.cellId);
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
            BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
                const CellSimilarityGraphVertex& vertex = graph[v];
                const string metaDataValue = getMetaData(vertex.cellId, metaDataName);
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
            BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
                CellSimilarityGraphVertex& vertex = graph[v];
                const string metaData = getMetaData(vertex.cellId, metaDataName);
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
                        colorString = brewerSetColor(iColor);
                    }
                    colorMap.insert(make_pair(group, colorString));
                }

            } else {

                // Each color gets used for a single meta data category.
                for(size_t group=0; group<sortedFrequencyTable.size(); group++) {
                    string colorString = "black";
                    if(group<12) {
                        colorString = brewerSetColor(group);
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
            BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
                CellSimilarityGraphVertex& vertex = graph[v];
                vertex.group = 0;
                vertex.color = getMetaData(vertex.cellId, metaDataName);
            }
        }



        // Color by meta data, interpreting the meta data value as a number.
        // We store in each vertex the meta data value that will determine the vertex color.
        else if(metaDataMeaning == "number") {
            colorByNumber = true;
            BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
                CellSimilarityGraphVertex& vertex = graph[v];
                vertex.group = 0;
                vertex.value = std::numeric_limits<double>::max();
                try {
                    vertex.value = lexical_cast<double>(getMetaData(vertex.cellId, metaDataName));
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
            BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
                CellSimilarityGraphVertex& vertex = graph[v];
                vertex.color.clear();
                vertex.group = 0;
            }
        }



        // Don't color the edges.
        BGL_FORALL_EDGES(e, graph, CellSimilarityGraph) {
            graph[e].color.clear();
        }

    }



    // Otherwise, all vertices and edges are colored black.
    else {

        BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
            CellSimilarityGraphVertex& vertex = graph[v];
            vertex.color.clear();
            vertex.group = 0;
        }
        BGL_FORALL_EDGES(e, graph, CellSimilarityGraph) {
            graph[e].color.clear();
        }
    }



    // If coloring by number, compute the color of each vertex.
    double minValue = std::numeric_limits<double>::max();
    double maxValue = std::numeric_limits<double>::lowest();
    if(colorByNumber) {

        // Compute the minimum and maximum values.
        BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
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
            BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
                graph[v].color = "black";
            }
        } else {
            const double scalingFactor = 1./(maxColorValue - minColorValue);
            BGL_FORALL_VERTICES(v, graph, CellSimilarityGraph) {
                CellSimilarityGraphVertex& vertex = graph[v];
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
    graph.writeSvg(html, hideEdges=="on", svgSizePixels, xViewBoxCenter, yViewBoxCenter, viewBoxHalfSize, vertexRadius, edgeThickness, colorMap);
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
            html << "<tr><td>Minimum<td class=centered>" << minValue;
            const int n = 4;
            for(int i=0; i<=n; i++) {
                const double x = i * (1./n);
                html << "<tr><td style='background-color:" << spectralColor(x) << "'>";
                html << "<td class=centered>";
                if(i==0) {
                    html << "<input type=text style='text-align:center;' form=coloringForm name=minColorValue id=minColorInput value='" << minColorValue << "'>";
                } else if(i==n) {
                    html << "<input type=text style='text-align:center;' form=coloringForm name=maxColorValue id=maxColorInput value='" << maxColorValue << "'>";
                } else {
                    html << minColorValue + x*(maxColorValue-minColorValue);;
                }
            }
            html << "<tr><td>Maximum<td class=centered>" << maxValue;
            html << "</table>";
        }
    }



    // End the div to contain the graphics and the table of meta data groups, side by side.
    html << "</div>";
}



void ExpressionMatrix::clusterDialog(
    const vector<string>& request,
    ostream& html)
{
    // Get the graph name.
    string graphName;
    if(!getParameterValue(request, "graphName", graphName)) {
        html << "Missing graph name.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }

    // Write the title.
    html << "<h1>Run clustering on graph " << graphName << "</h1>";

    // Form to enter the clustering parameters.
    html <<
        "<p>Enter parameters for label propagation clustering:"
        "<form action=cluster>"
        "Random number generator seed: <input type=text name=seed value=231>"
        "<br>Stop after this many iterations without changes: <input type=text name=stableIterationCountThreshold value=3>"
        "<br>Maximum number of iterations: <input type=text name=maxIterationCount value=100>"
		"<br>Meta data name to store the cluster of each cell (leave empty for none): <input type=text name=metaDataName>"
		"<br>Similarity threshold for graph edges: <input type=text name=similarityThreshold value='0.8'>"
        "<input type=hidden name=graphName value=" << graphName << ">"
        "<br><input type=submit value='Run clustering' autofocus>"
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
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
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
    double similarityThreshold;
    getParameterValue(request, "similarityThreshold", similarityThreshold);


    // Find the cell similarity graph.
    const auto it = graphs.find(graphName);
    if(it == graphs.end()) {
        html << "<p>Graph " << graphName << " does not exists.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }
    // const GraphCreationParameters& graphCreationParameters = it->second.first;
    CellSimilarityGraph& graph = *(it->second.second);

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
    ClusterGraph clusterGraph(graph);

	// Compute the average expression for each cluster - that is, for each vertex
    // of the cluster graph.
	BGL_FORALL_VERTICES(v, clusterGraph, ClusterGraph) {
		ClusterGraphVertex& vertex = clusterGraph[v];
		const size_t normalization = 2;	// Use L2 normalization. We might need to make this configurable.
		computeAverageExpression(vertex.cells, vertex.averageGeneExpression, normalization);
	}

	// Store in each edge the similarity of the two clusters, computed using the clusters
	// average expression stored in each vertex.
	clusterGraph.computeSimilarities();

	// Remove edges with low similarity.
	clusterGraph.removeWeakEdges(similarityThreshold);	// This may need to be made configurable.
    html << "<p>The cluster graph has " << num_vertices(clusterGraph);
    html << " vertices and " << num_edges(clusterGraph) << " edges.</p>";

	// Write out the cluster graph in graphviz format.
    clusterGraph.write("ClusterGraph.dot", geneNames);



    // Use Graphviz to create a layout of the ClusterGraph in svg format.
    // For better layouts, we would like to use -Goverlap=false because that does not work with the
    // ubuntu graphviz package (it is built without the triangulation library).
    const string command = "sfdp -O -T svg ClusterGraph.dot -Goverlap=scalexy -Gsplines=true";
    const int commandStatus = ::system(command.c_str());
    if(WIFEXITED(commandStatus)) {
    	const int exitStatus = WEXITSTATUS(commandStatus);
		if(exitStatus!=0 && exitStatus!=1) {	// sfdp returns 1 all the time just because of the message about missing triangulation.
			throw runtime_error("Error " + lexical_cast<string>(exitStatus) + " running graph layout command: " + command);
		}
    } else if(WIFSIGNALED(commandStatus)) {
    	const int signalNumber = WTERMSIG(commandStatus);
		throw runtime_error("Signal " + lexical_cast<string>(signalNumber) + " while running graph layout command: " + command);
    } else {
		throw runtime_error("Abnormal status " + lexical_cast<string>(commandStatus) + " while running graph layout command: " + command);

    }



    // Copy the output svg file to html.
    html << "<p>The cluster graph is shown below. Each vertex represents a cluster. "
    		"The most expressed genes for each cluster are listed."
    		"The size of each cluster is loosely related to the number of cells, but it is "
    		"also affected by the space needed to display the gene expression information.<br>";

    html << ifstream("ClusterGraph.dot.svg").rdbuf();

    // Add a button to continue.
    html << "<p><form action=graphs><input type=submit autofocus value=Continue></form>";
}



void ExpressionMatrix::createNewGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters.
    string graphName;
    if(!getParameterValue(request, "graphName", graphName)) {
        html << "Missing graph name.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }

    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "Missing cell set name.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }

    string similarPairsName;
    if(!getParameterValue(request, "similarPairsName", similarPairsName)) {
        html << "Missing similar pairs name.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }

    double similarityThreshold;
    if(!getParameterValue(request, "similarityThreshold", similarityThreshold)) {
        html << "Missing or invalid similarity threshold.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }

    int maxConnectivity;
    if(!getParameterValue(request, "maxConnectivity", maxConnectivity)) {
        html << "Missing or invalid max connectivity.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }

    // Check that the name does not already exist.
    if(graphs.find(graphName) != graphs.end()) {
        html << "<p>Graph " << graphName << " already exists.";
        html << "<p><form action=graphs><input type=submit value=Continue></form>";
        return;
    }


    // Create the graph.
    html << "<div style='font-family:courier'>";
    html << timestamp << "Graph creation begins.";
    createCellSimilarityGraph(graphName, cellSetName, similarPairsName, similarityThreshold, maxConnectivity);
    const GraphInformation& graphInfo = graphs[graphName].first;
    html <<
        "<br>" << timestamp << "New graph " << graphName << " was created. It has " << graphInfo.vertexCount <<
        " vertices and " << graphInfo.edgeCount << " edges"
		" after " << graphInfo.isolatedVertexCount << " isolated vertices were removed.";
    html << "</div>";

    // Add a button to continue.
    html << "<p><form action=graphs><input type=submit autofocus value=Continue></form>";


}



void ExpressionMatrix::removeGraph(
    const vector<string>& request,
    ostream& html)
{
    string graphName;
    if(!getParameterValue(request, "graphName", graphName)) {
        html << "<p>Missing graph name.";
    } else {
        const auto it = graphs.find(graphName);
        if(it == graphs.end()) {
            html << "<p>Graph " << graphName << " does not exist.";
        } else {
            graphs.erase(it);
            html << "<p>Graph " << graphName << " was removed.";
        }
    }


    html << "<p><form action=graphs><input type=submit autofocus value=Continue></form>";
}



// Get a list of the currently available sets of similar pairs.
void ExpressionMatrix::getAvailableSimilarPairs(
    vector<string>& availableSimilarPairs) const
{

    boost::regex regex(directoryName + "/SimilarPairs-(.*)-Info");
    boost::smatch matches;
    using boost::filesystem::directory_iterator;
    for(auto it=directory_iterator(directoryName); it!=directory_iterator(); ++it) {
        const string fileName = it->path().string();
        boost::smatch matches;
        if(!boost::regex_match(fileName, matches, regex)) {
            continue;
        }
        CZI_ASSERT(matches.size() == 2);
        availableSimilarPairs.push_back(matches[1]);
    }

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
        const string metaDataValue0 = getMetaData(cellId, metaDataNameId0);
        const string metaDataValue1 = getMetaData(cellId, metaDataNameId1);
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
    html << "<tr style='vertical-align:top;'><th>" << cellMetaDataNames[metaDataNameId1] << "&#129094";
    html << "<br>" << cellMetaDataNames[metaDataNameId0] << "&#129095";
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


