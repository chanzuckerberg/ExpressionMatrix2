// This file contains the implementation of http server functionality
// of class ExpressionMatrix related to cells.

#include "ExpressionMatrix.hpp"
#include "orderPairs.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "boost_tuple_tuple.hpp"



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
        html << " multiple title='Select two or more'";
    } else {
        html << " title='Select one'";
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



void ExpressionMatrix::exploreCell(
    const vector<string>& request,
    ostream& html)
{

    // Get the cell id.
    string cellIdString;
    const bool cellIdIsPresent = getParameterValue(request, "cellId", cellIdString);
    const CellId cellId = cellIdFromString(cellIdString);

    // Get the name of the gene set for which we want to write expression counts.
    string geneSetName = "AllGenes";
    getParameterValue(request, "geneSetName", geneSetName);

    // Locate the gene set.
    const auto it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        html << "<p>Gene set " << geneSetName << " does not exist.";
        return;
    }
    const auto& geneSet = it->second;

    // Write the form to get the cell id and the gene set.
    html <<
        "<form>"
        "<br><input type=text name=cellId autofocus>"
        " Specify a cell using a case-sensitive name or a numeric cell id between 0 and " << cellCount()-1 <<
        " included.<br>";
    writeGeneSetSelection(html, "geneSetName", false);
    html <<
        " Specify a gene set to display expression counts for this cell."
        "<br><input type=submit value=Display></form>";

    // If there is no cell id, do nothing.
    if(!cellIdIsPresent) {
        return;
    }

    // Access the cell.
    if(cellId==invalidCellId || cellId>cellCount()) {
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
    html << "<tr><td>Total number of genes with non-zero expression counts<td>" <<
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
    html << "<h2>Expression counts for this cell for gene set " << geneSetName << "</h2>";
    html <<
        "<p><strong>The following table of expression counts for this cell is sortable.</strong> Click on a header to sort by that header. "
        "Click again to reverse the sorting order."
        "<p><table id=countTable class=tablesorter><thead><tr><th>Gene<br>name<th>Raw<br>count"
        "<th>L1-normalized<br>count<br>(sum is 1)"
        "<th>L2-normalized<br>count<br>(sum<br>of<br>squares is 1)</thead><tbody>";



    if(geneSetName == "AllGenes") {

        // This is the old code that uses the global cell expression counts
        // stored in class ExpressionMatrix.
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
    } else {



        // This is the new code that uses computeExpressionVector to compute
        // the cell expression vector for this gene set.

        // Compute an unnormalized expression vector for this cell and gene set.
        vector< pair<GeneId, float> > expressionVector0;
        computeExpressionVector(cellId, geneSet, NormalizationMethod::none, expressionVector0);

        // Sort in decreasing order.
        sort(
            expressionVector0.begin(),
            expressionVector0.end(),
            OrderPairsBySecondGreaterThenByFirstLess< pair<GeneId, float> >());

        // Compute normalization factors.
        double sum1 = 0.;
        double sum2 = 0.;
        for(const auto& p: expressionVector0) {
            const double count = double(p.second);
            sum1 += count;
            sum2 += count*count;
        }
        const float factor1 = float(1./sum1);
        const float factor2 = float(1./sqrt(sum2));

        // Write it out.
        for(const auto& p: expressionVector0) {
            const GeneId localGeneId = p.first;
            const GeneId globalGeneId = geneSet.getGlobalGeneId(localGeneId);
            CZI_ASSERT(globalGeneId < geneCount());
            const string geneName = geneNames[globalGeneId];
            const float count = p.second;
            html <<  "<tr><td class=centered><a href=gene?geneId=" << urlEncode(geneName) << ">" << geneName << "</a>";
            html <<
                "<td class=centered>" << count;
            const auto oldPrecision = html.precision(3);
            html <<
                "<td class=centered>" << count * factor1 <<
                "<td class=centered>" << count * factor2;
            html.precision(oldPrecision);
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



void ExpressionMatrix::compareTwoCells(
    const vector<string>& request,
    ostream& html)
{
    // Get the cell ids.
    array<string, 2> cellIdStrings;
    const bool cellId0IsPresent = getParameterValue(request, "cellId0", cellIdStrings[0]);
    const bool cellId1IsPresent = getParameterValue(request, "cellId1", cellIdStrings[1]);
    const bool cellIdsArePresent = cellId0IsPresent && cellId1IsPresent;

    // Get the name of the gene set to be used for the comparison.
    string geneSetName = "AllGenes";
    getParameterValue(request, "geneSetName", geneSetName);


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
    html << "><br>Select a gene set to be used for the comparison: ";
    writeGeneSetSelection(html, "geneSetName", {geneSetName}, false);
    html <<
        "<br><input type=submit value=Compare>"
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

    // Access the gene set.
    const auto it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        html << "<p>Gene set " << geneSetName << " does not exist." << endl;
        return;
    }
    const GeneSet& geneSet = it->second;

    // Write a title.
    html  << "<h1>Comparison of cells " << cellIds[0] << " ";
    writeCellLink(html, cellIds[0], false);
    html << " and " << cellIds[1] << " ";
    writeCellLink(html, cellIds[1], false);
    html << "using gene set " << geneSetName << "</h1>";


    // Write a table of similarities between these two cells.
    html <<
        "<p>The similarity of these two cells "
        "computed using gene set " << geneSetName << " is ";
    const auto oldPrecision = html.precision(3);
    html << computeCellSimilarity(geneSet, cellIds[0], cellIds[1]);
    html.precision(oldPrecision);
    html << ".<br>";



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
                const GeneId geneId = it1->first;
                if(geneSet.contains(geneId)) {
                    const double count1 = it1->second;
                    data.push_back(make_tuple(count1, geneId, 0, count1));
                }
                ++it1;
            }
        } else {
            if(it1==end1) {
                // it1 is at end, it0 is not.
                const GeneId geneId = it0->first;
                if(geneSet.contains(geneId)) {
                    const double count0 = it0->second;
                    data.push_back(make_tuple(count0, geneId, count0, 0));
                }
                ++it0;
            } else {
                // Neither it0 nor it1 are at end.
                if(it0->first < it1->first) {
                    const GeneId geneId = it0->first;
                    if(geneSet.contains(geneId)) {
                        const double count0 = it0->second;
                        data.push_back(make_tuple(count0, geneId, count0, 0));
                    }
                    ++it0;
                } else if(it1->first < it0->first) {
                    const GeneId geneId = it1->first;
                    if(geneSet.contains(geneId)) {
                        const double count1 = it1->second;
                        data.push_back(make_tuple(count1, geneId, 0, count1));
                    }
                    ++it1;
                } else {
                    const GeneId geneId = it0->first;
                    if(geneSet.contains(geneId)) {
                        const double count0 = it0->second;
                        const double count1 = it1->second;
                        data.push_back(make_tuple(count0+count1, geneId, count0, count1));
                    }
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
        "<input type=button value='Zoom in' onclick='scale(0.5);' style='width:20em;border-radius:5px;' /><br>"
        "<input type=button value='Zoom out' onclick='scale(2);' style='width:20em;border-radius:5px;' /><br>"
        "<input type=button value='Horizontal only zoom in' onclick='horizontalScale(0.5);' style='width:20em;border-radius:5px;' /><br>"
        "<input type=button value='Horizontal only zoom out' onclick='horizontalScale(2);' style='width:20em;border-radius:5px;' /><br>"
        "<input type=button value='Vertical only zoom in' onclick='verticalScale(0.5);' style='width:20em;border-radius:5px;' /><br>"
        "<input type=button value='Vertical only zoom out' onclick='verticalScale(2);' style='width:20em;border-radius:5px;' /><br>"
        "<input type=button value='Larger marker' onclick='changeMarkerSize(1);' style='width:20em;border-radius:5px;' /><br>"
        "<input type=button value='Smaller marker' onclick='changeMarkerSize(-1);' style='width:20em;border-radius:5px;' /><br>"
        "</div>"
        "<div style='clear:both;' />"
        ;



    // Write to html jQuery and TableSorter so we can make the table below sortable.
    writeJQuery( html);
    writeTableSorter(html);



    // Write the table with the counts.
    html <<
        "<p><strong>The table is sortable.</strong> Click on a header to sort by that header. "
        "Click again to reverse the sorting order."
        "<table id=expressionTable class=tablesorter><thead>"
        "<tr><th>Gene<br>id<th>Gene<br>name<th>Count for<br>cell " << cellIds[0] << "<br>";
    writeCellLink(html, cellIds[0], false);
    html << "<th>Count for<br>cell " << cellIds[1] << "<br>";
    writeCellLink(html, cellIds[1], false);
    html << "</thead><tbody>";

    for(const auto& t: data) {
        const int geneId = t.get<1>();
        const double count0 = t.get<2>();
        const double count1 = t.get<3>();
        html << "<tr><td>";
        writeGeneLink(html, geneId, true);
        html << "<td>";
        writeGeneLink(html, geneId, false);
        html << "<td class=centered>" << count0 << "<td class=centered>" << count1;
    }
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#expressionTable').tablesorter();"
        "window.location='#clusterGraph'});"
        "</script>"
        ;


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


    // Form to create a new cell set from string meta data.
    html <<
        "<br><h2>Create a new cell set using string meta data</h2>"
        "<p><form action=createCellSetUsingMetaData>"
        "<input type=submit value='Create a new cell set'> named "
        "<input type=text required name=cellSetName>"
        " consisting of cells for which meta data field ";
        writeMetaDataSelection(html, "metaData", false);
    html <<
        " matches this "
        "<select name=matchType>"
        "<option value=stringMatch selected=selected>string</option>"
        "<option value=regexMatch>regular expression</option>"
        "</select>:"
        "<input type=text name=matchTarget>"
        "</form>";



    // Form to create a new cell set from numeric meta data.
    html <<
        "<br><h2>Create a new cell set using numeric meta data</h2>"
        "<p><form action=createCellSetUsingNumericMetaData>"
        "<input type=submit value='Create a new cell set'> named "
        "<input type=text required name=cellSetName>"
        " consisting of cells for which meta data field ";
        writeMetaDataSelection(html, "metaData", false);
    html <<
        " is numeric, greater than "
        "<input type=text name=lowerBound>"
        " (leave blank to ignore), and less than "
        "<input type=text name=upperBound> (leave blank to ignore)"
        "</form>";



    // Form to create a new cell set by union/intersection existing cell sets.
    html <<
        "<br><h2>Create a new cell set by union/intersection of existing cell sets</h2>"
        "<p><form action=createCellSetIntersectionOrUnion>"
        "<input type=submit value='Create a new cell set'> named "
        "<input type=text required name=cellSetName>"
        " as the "
        "<select name=operation>"
        "<option value=union>union</option>"
        "<option value=intersection>intersection</option>"
        "</select>"
        " of the selected cell sets: ";
    writeCellSetSelection(html, "inputCellSets", true);
    html << "</form>";



    // Form to create a new cell set as the set difference of existing cell sets.
    html <<
        "<br><h2>Create a new cell set as the set difference of existing cell sets</h2>"
        "<p><form action=createCellSetDifference>"
        "<input type=submit value='Create a new cell set'> named "
        "<input type=text required name=cellSetName>"
        " as the set difference of cell set ";
    writeCellSetSelection(html, "inputCellSet0", false);
    html << " minus cell set ";
    writeCellSetSelection(html, "inputCellSet1", false);
    html << ".</form>";



    // Form to create a new cell set by downsampling an existing cell set.
    html <<
        "<br><h2>Create a new cell set by downsampling an existing cell set</h2>"
        "<p>The new cell set will be a random subset of the specified cell set."
        " Each cell in the specified cell set is inserted in the random subset with the specified probability."
        " Therefore, the downsampling rate will be approximately equal to the specified probability."
        "<p><form action=downsampleCellSet>"
        "<input type=submit value='Create a new cell set'> named "
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



// Get a copy of the cell set with a given name, for use in the Python API.
// Returns an empty cell set if a cel set with the specified name does not exist.
vector<CellId> ExpressionMatrix::getCellSet(const string& cellSetName) const
{
    vector<CellId> cells;
    const auto it = cellSets.cellSets.find(cellSetName);
    if(it != cellSets.cellSets.end()) {
        const CellSet& cellSet = *(it->second);
        cells.resize(cellSet.size());
        copy(cellSet.begin(), cellSet.end(), cells.begin());
    }
    return cells;
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

    // Get the string or reggular expression to be matches.
    string matchTarget;
    if(!getParameterValue(request, "matchTarget", matchTarget)) {
        html << "Missing string or regular expression to match.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }
    string decodedMatchTarget;
    urlDecode(matchTarget, decodedMatchTarget);

    // Figure out if we are supposed to match a string or a regular expression.
    string matchType = "stringMatch";
    getParameterValue(request, "matchType", matchType);
    const bool useRegex = (matchType == "regexMatch");

    if(createCellSetUsingMetaData(cellSetName, metaData, decodedMatchTarget, useRegex)) {
        html << "<p>Newly created cell set " << cellSetName << " has ";
        html << cellSets.cellSets[cellSetName]->size() << " cells.";
    } else {
        html << "<p>Unable to create cell set " << cellSetName << ".";
    }
    html << "<p><form action=cellSets><input type=submit value=Continue></form>";

}



void ExpressionMatrix::createCellSetUsingNumericMetaData(const vector<string>& request, ostream& html)
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

    bool useLowerBound = false;
    double lowerBound = 0.;
    string lowerBoundString;
    if(getParameterValue(request, "lowerBound", lowerBoundString)) {
        try {
            lowerBound = lexical_cast<double>(lowerBoundString);
            useLowerBound = true;
        } catch(bad_lexical_cast) {
            html << "Invalid lower bound value.";
        }
    }

    bool useUpperBound = false;
    double upperBound = 0.;
    string upperBoundString;
    if(getParameterValue(request, "upperBound", upperBoundString)) {
        try {
            upperBound = lexical_cast<double>(upperBoundString);
            useUpperBound = true;
        } catch(bad_lexical_cast) {
            html << "Invalid upper bound value.";
        }
    }


    try {
        createCellSetUsingNumericMetaData(
            cellSetName, metaData,
            useLowerBound, lowerBound,
            useUpperBound, upperBound);
        html << "<p>Newly created cell set " << cellSetName << " has ";
        html << cellSets.cellSets[cellSetName]->size() << " cells.";
    } catch(...) {
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



void ExpressionMatrix::createCellSetDifference(const vector<string>& request, ostream& html)
{
    // Get the name of the cell set to be created.
    string cellSetName;
    if(!getParameterValue(request, "cellSetName", cellSetName)) {
        html << "Missing cell set name.";
        html << "<p><form action=cellSets><input type=submit value=Continue></form>";
        return;
    }



    // Get the names of the input cell sets.
    string inputCellSet0, inputCellSet1;
    getParameterValue(request, "inputCellSet0", inputCellSet0);
    getParameterValue(request, "inputCellSet1", inputCellSet1);



    // Do the difference.
    if(createCellSetDifference(inputCellSet0, inputCellSet1, cellSetName)) {
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
    double probability = 0.1;
    getParameterValue(request, "probability", probability);
    int seed = 231;
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
    try {
    	filesystem::remove(fileName);
        html << "<p>Cell set " << cellSetName << " was removed.";
    } catch(...) {
        html << "<p>Cell set " << cellSetName << " was removed from memory but the corresponding memory mapped file ";
        html << fileName << " could not be removed.";
    }
    html << "<p><form action=cellSets><input type=submit value=Continue></form>";

}
