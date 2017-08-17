// This file contains the implementation of http server functionality
// of class ExpressionMatrix related to genes.

#include "ExpressionMatrix.hpp"
#include "tokenize.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



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



ostream& ExpressionMatrix::writeGeneSetSelection(
    ostream& html,
    const string& selectName,
    bool multiple) const
{
    return writeGeneSetSelection(html, selectName, {"AllGenes"}, multiple);
}



ostream& ExpressionMatrix::writeGeneSetSelection(
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
    for(const auto& p: geneSets) {
        const string& geneSetName = p.first;
        html << "<option value=" << geneSetName;
        if(selected.find(geneSetName) != selected.end()) {
            html << " selected=selected";
        }
        html << ">" << geneSetName << "</option>";
    }
    html << "</select>";

    return html;
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

    // Write the form for information on a single gene.
    html <<
        "<form>"
        "<input type=submit value='Show information about gene'>"
        " <input type=text name=geneId";
    if(geneIdIsPresent) {
        html << " value=" << geneIdString;
    }
    html <<
        " autofocus>"
        " (specify a gene using a case-sensitive name or a numeric id between 0 and  " << geneCount()-1 << "),";
    html << "<br>displaying expression counts for this gene for cells in cell set ";
    const set<string> selectedCellSet = {cellSetName};
    writeCellSetSelection(html, "cellSetName", selectedCellSet, false);
    html << "<br>and showing for each cell the following cell meta data fields: ";
    writeMetaDataSelection(html, "metaDataName", metaDataToDisplay, true);
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
            html << "<td class=centered>" << getCellMetaData(data.cellId, metaDataNameStringId);
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



void ExpressionMatrix::createGeneSetFromRegex(const vector<string>& request, ostream& html)
{
    string geneSetName;
    if(!getParameterValue(request, "geneSetName", geneSetName)) {
        html << "Missing gene set name.";
        html << "<p><form action=geneSets><input type=submit value=Continue></form>";
        return;
    }

    string regex;
    if(!getParameterValue(request, "regex", regex)) {
        html << "Missing regular expression.";
        html << "<p><form action=geneSets><input type=submit value=Continue></form>";
        return;
    }
    string decodedRegex;
    urlDecode(regex, decodedRegex);

    if(createGeneSetFromRegex(geneSetName, decodedRegex)) {
        html << "<p>Newly created gene set " << geneSetName << " has ";
        html << geneSets[geneSetName].size() << " genes.";
    } else {
        html << "<p>Unable to create gene set " << geneSetName << ".";
    }
    html << "<p><form action=geneSets><input type=submit value=Continue></form>";
}



void ExpressionMatrix::createGeneSetFromGeneNames(const vector<string>& request, ostream& html)
{
    string geneSetName;
    if(!getParameterValue(request, "geneSetName", geneSetName)) {
        html << "Missing gene set name.";
        html << "<p><form action=geneSets><input type=submit value=Continue></form>";
        return;
    }

    string geneNamesString;
    if(!getParameterValue(request, "genes", geneNamesString)) {
        html << "Missing gene names.";
        html << "<p><form action=geneSets><input type=submit value=Continue></form>";
        return;
    }

    // Parse the gene names into a vector.
    vector<string> geneNames;
    tokenize(" ,\t\n\r", geneNamesString, geneNames, false);

#if 0
    // Debug output to show detail of tokenization.
    html << "<p>";
    for(const char c: geneNamesString) {
        html << hex << int(c) << dec << " ";
    }
    for(const string& geneName: geneNames) {
        html << "<p>Gena name of length " << geneName.size() << ": ";
        for(const char c: geneName) {
            html << hex << int(c) << dec << " ";
        }
    }
#endif

    // Create the gene set.
    int ignoredCount, emptyCount;
    if(createGeneSetFromGeneNames(geneSetName, geneNames, ignoredCount, emptyCount)) {
        html << "<p>Newly created gene set " << geneSetName << " has ";
        html << geneSets[geneSetName].size() << " genes.";
        if(ignoredCount) {
            html << "<p>" << ignoredCount << " of the " << geneNames.size()-emptyCount <<
                " specified names were ignored because they don't correspond to valid gene names.";
        }
    } else {
        html << "<p>Unable to create gene set " << geneSetName << ".";
    }
    html << "<p><form action=geneSets><input type=submit value=Continue></form>";
}



void ExpressionMatrix::createGeneSetIntersectionOrUnion(const vector<string>& request, ostream& html)
{
    // Get the name of the gene set to be created.
    string geneSetName;
    if(!getParameterValue(request, "geneSetName", geneSetName)) {
        html << "Missing gene set name.";
        html << "<p><form action=geneSets><input type=submit value=Continue></form>";
        return;
    }

    // Get the name of the operation to be performed (intersection or union).
    string operation;
    if(!getParameterValue(request, "operation", operation)) {
        html << "Missing operation.";
        html << "<p><form action=geneSets><input type=submit value=Continue></form>";
        return;
    }
    bool doUnion;
    if(operation == "intersection") {
        doUnion = false;
    } else if(operation == "union") {
        doUnion = true;
    } else {
        html << "Invalid operation.";
        html << "<p><form action=geneSets><input type=submit value=Continue></form>";
        return;
    }


    // Get the names of the input gene sets.
    set<string> inputGeneSets;
    getParameterValues(request, string("inputGeneSets"), inputGeneSets);
    if(inputGeneSets.size() < 2) {
        html << "At least two input gene sets should be specified.";
        return;
    }

    // Concatenate the input gene sets with commas.
    string inputGeneSetsString;
    for(const string& inputGeneSet: inputGeneSets) {
        inputGeneSetsString.append(inputGeneSet);
        inputGeneSetsString.append(",");
    }
    inputGeneSetsString.resize(inputGeneSetsString.size()-1);


    // Do the intersection or union.
    if(createGeneSetIntersectionOrUnion(inputGeneSetsString, geneSetName, doUnion)) {
        html << "<p>Newly created gene set " << geneSetName << " has ";
        html << geneSets[geneSetName].size() << " genes.";
    } else {
        html << "<p>Unable to create gene set " << geneSetName << ".";
    }
    html << "<p><form action=geneSets><input type=submit value=Continue></form>";

}



void ExpressionMatrix::createGeneSetDifference(const vector<string>& request, ostream& html)
{
    // Get the name of the gene set to be created.
    string geneSetName;
    if(!getParameterValue(request, "geneSetName", geneSetName)) {
        html << "Missing gene set name.";
        html << "<p><form action=geneSets><input type=submit value=Continue></form>";
        return;
    }



    // Get the names of the input gene sets.
    string inputGeneSet0, inputGeneSet1;
    getParameterValue(request, "inputGeneSet0", inputGeneSet0);
    getParameterValue(request, "inputGeneSet1", inputGeneSet1);



    // Do the difference.
    if(createGeneSetDifference(inputGeneSet0, inputGeneSet1, geneSetName)) {
        html << "<p>Newly created gene set " << geneSetName << " has ";
        html << geneSets[geneSetName].size() << " genes.";
    } else {
        html << "<p>Unable to create gene set " << geneSetName << ".";
    }
    html << "<p><form action=geneSets><input type=submit value=Continue></form>";
}



void ExpressionMatrix::createGeneSetUsingInformationContent(const vector<string>& request, ostream& html)
{
    // Get the name of the gene set to use to compute gene information content.
    string geneSetName;
    getParameterValue(request, "geneSetName", geneSetName);

    // Get the name of the new gene set to be created.
    string newGeneSetName;
    getParameterValue(request, "newGeneSetName", newGeneSetName);

    // Get the name of the cell set to use to compute gene information content.
    string cellSetName;
    getParameterValue(request, "cellSetName", cellSetName);

    // Get the normalization method to be used to compute gene information content.
    const NormalizationMethod normalizationMethod = getNormalizationMethod(request, NormalizationMethod::L2);

    // Get the threshold on gene information content.
    // This is the minimum gene information content required for a gene to be included
    // in the gene set being created.
    float threshold = 2;
    getParameterValue(request, "threshold", threshold);

    // Locate the gene set.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        html << "<p>Invalid gene set name " << geneSetName;
        return;
    }
    const GeneSet& geneSet = itGeneSet->second;

    // Verify that the new gene set does not already exist.
    if(geneSets.find(newGeneSetName) != geneSets.end()) {
        html << "<p>Gene set " << newGeneSetName << " already exists.";
        return;
    }

    // Locate the cell set.
    const auto itCellSet = cellSets.cellSets.find(cellSetName);
    if(itCellSet == cellSets.cellSets.end()) {
        html << "<p>Invalid cell set name " << cellSetName;
        return;
    }
    const CellSet& cellSet = *(itCellSet->second);

    // Write a title.
    html << "<h1>Creation of gene set " << newGeneSetName;
    html << " using gene information content</h1>";

    // Write a table summarizing the creation parameters used.
    html <<
        "<p>Summary of creation parameters for new gene set " << newGeneSetName <<
        ":<p><table style='table-layout:fixed;width:500px;'>"
        "<tr><th class=left style='width:300px;'>Input gene set<td class=centered>" << geneSetName <<
        " (" << geneSet.size() << " genes)"
        "<tr><th class=left style='width:300px;'>Cell set used to compute gene information content<td class=centered>" <<
        cellSetName << " (" << cellSet.size() << " cells)"
        "<tr><th class=left style='width:300px;'>Normalization method used to compute gene information content<td class=centered>" <<
        normalizationMethodToLongString(normalizationMethod) <<
        "<tr><th class=left style='width:300px;'>Gene information content threshold (bits)<td class=centered>" << threshold <<
        "</table>"
        "<p>This preliminary implementation is very slow (typically around 10 seconds per thousand cells)."
        " Future versions will be faster.</p>";
        ;

    // Compute gene information content using the requested normalization method.
    vector<float> informationContent;
    computeGeneInformationContent(geneSet, cellSet, normalizationMethod, informationContent);

    // Create the new gene set.
    GeneSet& newGeneSet = geneSets[newGeneSetName];
    newGeneSet.createNew(directoryName + "/GeneSet-" + newGeneSetName);
    for(GeneId localGeneId=0; localGeneId!=geneSet.size(); localGeneId++) {
        if(informationContent[localGeneId] > threshold) {
            const GeneId globalGeneId = geneSet.getGlobalGeneId(localGeneId);
            newGeneSet.addGene(globalGeneId);
        }
    }
    html << "<p>Gene set " << newGeneSetName << " created. It has " << newGeneSet.size() << " genes.";


    html << "<p><form action=geneSets><input type=submit value=Continue></form>";
}



void ExpressionMatrix::exploreGeneSets(
    const vector<string>& request,
    ostream& html)
{
    // Write a title.
    html << "<h1>Gene sets</h1>";

    // Write a table listing the gene sets in existence.
    html << "<p><table><th>Gene<br>set<br>name<th>Number<br>of<br>genes<th class=centered>Click<br>to<br>remove";
    for(const auto& p: geneSets) {
        const string& name = p.first;
        const auto& geneSet = p.second;
        html << "<tr><td><a href='geneSet?geneSetName=" << urlEncode(name) << "'>" << name << "</a><td class=centered>" << geneSet.size();
        html << "<td  class=centered>";
        if(name != "AllGenes") {
            html << "<a href='removeGeneSet?geneSetName=" << urlEncode(name) << "'>Remove</a>";
        }
    }
    html << "</table>";



    // Form to create a gene set from a regular expression.
    html <<
        "<br><h2>Create a new gene set using a regular expression for gene names</h2>"
        "<p><form action=createGeneSetFromRegex>"
        "<input type=submit value='Create a new gene set'> named "
        "<input type=text required name=geneSetName>"
        " consisting of genes with names matching this regular expression: "
        "<input type=text name=regex>"
        "</form>";



    // Form to create a gene set by pasting gene names.
    html <<
        "<br><h2>Create a new gene by pasting gene names</h2>"
        "<p><form id=createGeneSetFromGeneNames action=createGeneSetFromGeneNames>"
        "<input type=submit value='Create a new gene set'> named "
        "<input type=text required name=geneSetName>"
        " consisting of the following genes: "
        "</form>"
        "<textarea form=createGeneSetFromGeneNames name=genes rows=6 cols=80></textarea>"
        "<br>Paste gene names here. White space and commas can be used as separators.";



    // Form to create a new gene set by union/intersection of existing gene sets.
    html <<
        "<br><h2>Create a new gene set by union/intersection of existing gene sets</h2>"
        "<p><form action=createGeneSetIntersectionOrUnion>"
        "<input type=submit value='Create a new gene set'> named "
        "<input type=text required name=geneSetName>"
        " as the "
        "<select name=operation>"
        "<option value=union>union</option>"
        "<option value=intersection>intersection</option>"
        "</select>"
        " of the selected gene sets: ";
    writeGeneSetSelection(html, "inputGeneSets", true);
    html << "</form>";



    // Form to create a new gene set as the set difference of existing gene sets.
    html <<
        "<br><h2>Create a new gene set as the set difference of existing gene sets</h2>"
        "<p><form action=createGeneSetDifference>"
        "<input type=submit value='Create a new gene set'> named "
        "<input type=text required name=geneSetName>"
        " as the set difference of gene set ";
    writeGeneSetSelection(html, "inputGeneSet0", false);
    html << " minus gene set ";
    writeGeneSetSelection(html, "inputGeneSet1", false);
    html << ".</form>";



    // Form to create a gene set using gene information content.
    html <<
        "<br><h2>Create a gene set using gene information content</h2>"
        "<form action=createGeneSetUsingInformationContent>"
        "<input type=submit value='Create a new gene set using gene information content'>"
        " as follows: <br>starting with with gene set ";
    writeGeneSetSelection(html, "geneSetName", false);
    html << ",<br>compute gene information content using cell set ";
    writeCellSetSelection(html, "cellSetName", {"AllCells"}, false);
    html << " and ";
    writeNormalizationSelection(html, NormalizationMethod::L2);
    html <<
        ",<br>then keep genes with information content at least "
        "<input type=text name=threshold>"
        " bits<br>and name the resulting gene set "
        "<input type=text name=newGeneSetName>"
        "</form>";



    // Write the form to display gene information content.
    html <<
        "<h2>Show gene information content</h2>"
        "<form action=geneInformationContent>"
        "<input type=submit value='Show gene information content for gene set'> ";
    writeGeneSetSelection(html, "geneSetName", false);
    html << " cell set ";
    writeCellSetSelection(html, "cellSetName", {"AllCells"}, false);
    html << "</form>";

}



void ExpressionMatrix::exploreGeneSet(
    const vector<string>& request,
    ostream& html)
{
    // Get the name of the gene set we want to look at.
    string geneSetName;
    if(!getParameterValue(request, "geneSetName", geneSetName)) {
        html << "Missing gene set name.";
        return;
    }

    // Write a title.
    html << "<h1>Gene set " << geneSetName << "</h1>";

    // Locate the gene set.
    const auto it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        html << "<p>Gene set " << geneSetName << " does not exist.";
        return;
    }
    const auto& geneSet = it->second;
    html << "<p>This gene set has " << geneSet.size() << " genes." << endl;

    // Write a table containing the genes in this gene set.
    html << "<table>";
    for(const GeneId geneId: geneSet) {
        const string geneName = geneNames[geneId];
        html <<  "<tr><td class=centered><a href=gene?geneId=" << urlEncode(geneName) << ">" << geneName << "</a>";
    }
    html << "</table>";

}



void ExpressionMatrix::removeGeneSet(
    const vector<string>& request,
    ostream& html)
{
    // Get the name of the gene set we want to remove.
    string geneSetName;
    if(!getParameterValue(request, "geneSetName", geneSetName)) {
        html << "Missing gene set name.";
        return;
    }

    // Locate the gene set.
    const auto it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        html << "<p>Gene set " << geneSetName << " does not exist.";
        return;
    }

    // Remove it.
    GeneSet& geneSet = it->second;
    try{
        geneSet.remove();
    } catch(...) {
        html << "<p>Unable to remove gene set " << geneSetName << ".";
        throw;
    }
    geneSets.erase(it);
    html << "<p>Gene set " << geneSetName << " was removed.";

    html << "<p><form action=geneSets><input type=submit value=Continue></form>";
}



void ExpressionMatrix::exploreGeneInformationContent(const vector<string>& request, ostream& html)
{
    // Get the name of the gene set to use to compute gene information content.
    string geneSetName;
    getParameterValue(request, "geneSetName", geneSetName);

    // Get the name of the cell set to use to compute gene information content.
    string cellSetName;
    getParameterValue(request, "cellSetName", cellSetName);

    // Locate the gene set.
    const auto itGeneSet = geneSets.find(geneSetName);
    if(itGeneSet == geneSets.end()) {
        html << "<p>Invalid gene set name " << geneSetName;
        return;
    }
    const GeneSet& geneSet = itGeneSet->second;

    // Locate the cell set.
    const auto itCellSet = cellSets.cellSets.find(cellSetName);
    if(itCellSet == cellSets.cellSets.end()) {
        html << "<p>Invalid cell set name " << cellSetName;
        return;
    }
    const CellSet& cellSet = *(itCellSet->second);


    html <<
        "<h1>Gene information content</h1>"
        "<p>Gene information content computed taking only into account genes in gene set " << geneSetName <<
        " (" << geneSet.size() << " genes) "
        " and cells in cell set " << cellSetName <<
        " (" << cellSet.size() << " cells)."
        "<p><strong>The table below is sortable.</strong> Click on a header to sort by that column. "
        "Click again to reverse the sorting order."
        "<p>This preliminary implementation is very slow (typically around 20 seconds per thousand cells)."
        " Future versions will be faster.";

    // Compute gene information content for all normalization methods.
    vector<float> informationContent0;
    vector<float> informationContent1;
    vector<float> informationContent2;
    computeGeneInformationContent(geneSet, cellSet, NormalizationMethod::None, informationContent0);
    computeGeneInformationContent(geneSet, cellSet, NormalizationMethod::L1  , informationContent1);
    computeGeneInformationContent(geneSet, cellSet, NormalizationMethod::L2  , informationContent2);

    // Write to html jQuery and TableSorter so we can make the table below sortable.
    writeJQuery( html);
    writeTableSorter(html);

    // Write the table of gene information content.
    const auto oldPrecision = html.precision(3);
    html << "<table id=countTable class=tablesorter style='table-layout:fixed;width:480px;'><thead><th>Gene";
    for(const NormalizationMethod normalizationMethod: validNormalizationMethods) {
        html << "<th>Gene information content in bits computed using " << normalizationMethodToLongString(normalizationMethod);
    }
    html << "</thead><tbody>";
    for(GeneId localGeneId=0; localGeneId<geneSet.size(); localGeneId++) {
        const GeneId globalGeneId = geneSet.getGlobalGeneId(localGeneId);
        CZI_ASSERT(globalGeneId < geneCount());
        const string geneName = geneNames[globalGeneId];
        html <<  "<tr><td class=centered style='width:160px;'><a href='gene?geneId=" << urlEncode(geneName) << "'>" << geneName << "</a>";
        html <<
            "<td class=centered style='width:160px;'>" << informationContent0[localGeneId] <<
            "<td class=centered style='width:160px;'>" << informationContent1[localGeneId] <<
            "<td class=centered style='width:160px;'>" << informationContent2[localGeneId];
    }
    html.precision(oldPrecision);


    // Finish the table and make it sortable.
    html <<
        "</tbody></table>"
        "<script>"
        "$(document).ready(function(){$('#countTable').tablesorter({sortList:[[3,0]]});});"
        "</script>"
        ;
}
