// Http server functionality related to signature graphs.

#include "ExpressionMatrix.hpp"
#include "SignatureGraph.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



void ExpressionMatrix::exploreSignatureGraphs(const vector<string>& request, ostream& html)
{
    html << "<h1>Signature graphs</h1>";
    html << "<p>In a signature graph, all cells with the same LSH signature "
        "are collapsed into a single vertex.";

    // Table of existing signature graphs.
    html << "<h2>Existing signature graphs</h2><table>";
    for(const auto& p: signatureGraphs) {
        html << "<tr><td><a href='exploreSignatureGraph?signatureGraphName="
            << p.first << "'>" << p.first << "</a>"
            "<td><form action=removeSignatureGraph>"
            "<input type=submit value='Remove'>"
            "<input hidden type=text name=signatureGraphName value='" << p.first << "'>"
            "</form>";
    }
    html << "</table>";

    // Form to create a new signature graph.
    html <<
        "<h2>Create a new signature graph</h2>"
        "<form action=createSignatureGraph>"
        "<input type=submit value='Create'> a new signature graph named "
        "<input type=text size=8 name=signatureGraphName required> "
        "for cell set ";
    writeCellSetSelection(html, "cellSetName", false);
    html << " using LSH signatures ";
    writeLshSelection(html, "lshName");
    html << ". Only generate vertices for signatures with at least "
        "<input type=text size=8 name=minCellCount value=3 required> cells."
        "</form>";
}



void ExpressionMatrix::exploreSignatureGraph(const vector<string>& request, ostream& html)
{
    string signatureGraphName;
    getParameterValue(request, "signatureGraphName", signatureGraphName);
    if(signatureGraphName.empty()) {
        html << "Signature graph name is missing.";
        return;
    }
    SignatureGraph& signatureGraph = getSignatureGraph(signatureGraphName);

    html << "<h1>Signature graph " << signatureGraphName << "</h1>";

    SignatureGraph::SvgParameters svgParameters = signatureGraph.getDefaultSvgParameters();
    signatureGraph.writeSvg(html, svgParameters);

}



void ExpressionMatrix::createSignatureGraph(const vector<string>& request, ostream& html)
{
    string signatureGraphName;
    getParameterValue(request, "signatureGraphName", signatureGraphName);
    if(signatureGraphName.empty()) {
        html << "Signature graph name is missing.";
        return;
    }

    string lshName;
    getParameterValue(request, "lshName", lshName);
    if(lshName.empty()) {
        html << "Lsh name is missing.";
        return;
    }

    string cellSetName;
    getParameterValue(request, "cellSetName", cellSetName);
    if(cellSetName.empty()) {
        html << "Cell set name is missing.";
        return;
    }

    CellId minCellCount = 3;
    getParameterValue(request, "minCellCount", minCellCount);

    html << "<h1>Create signature graph " << signatureGraphName << "</h1>";
    createSignatureGraph(signatureGraphName, cellSetName, lshName, minCellCount);
    html << "<p>Signature graph " << signatureGraphName << " was created."
        "<p><form action=exploreSignatureGraph>"
        "<input type=text hidden name=signatureGraphName value=" << signatureGraphName <<
        "><input type=submit value=Continue></form>";

}



void ExpressionMatrix::removeSignatureGraph(const vector<string>& request, ostream& html)
{
    string signatureGraphName;
    getParameterValue(request, "signatureGraphName", signatureGraphName);
    if(signatureGraphName.empty()) {
        html << "Signature graph name is missing.";
        return;
    }
    removeSignatureGraph(signatureGraphName);

    html << "<p>Signature graph " << signatureGraphName << " was removed."
        "<p><form action=exploreSignatureGraphs>"
        "<input type=submit value=Continue></form>";

}
