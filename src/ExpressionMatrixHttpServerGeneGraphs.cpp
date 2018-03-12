// Http server functionality related to gene graphs.

#include "ExpressionMatrix.hpp"
#include "GeneGraph.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



void ExpressionMatrix::exploreGeneGraphs(const vector<string>& request, ostream& html)
{
    html << "<h1>Gene graphs</h1>";

    // Table of existing gene graphs.
    html << "<h2>Existing gene graphs</h2>";
    if(geneGraphs.empty()) {
        html << "<p>No gene graphs exist. You can create one using the form below.";
    } else {
        html << "<table>";
        for(const auto& p: geneGraphs) {
            html << "<tr><td><a href='exploreGeneGraph?geneGraphName="
                << p.first << "'>" << p.first << "</a>"
                "<td><form action=removeGeneGraph>"
                "<input type=submit value='Remove'>"
                "<input hidden type=text name=geneGraphName value='" << p.first << "'>"
                "</form>";
        }
        html << "</table>";
    }



    // Form to create a new gene graph.
    html <<
        "<h2>Create a new gene graph</h2>"
        "<form action=createGeneGraph>"
        "<table>"
        "<tr><td>Gene set<td class=centered>";
    writeGeneSetSelection(html, "geneSetName", false);
    html << "<tr><td>Similar gene pairs<td class=centered>";
    writeSimilarGenePairsSelection(html, "similarGenePairsName");
    html <<
        "<tr><td>Similarity threshold<td class=centered><input type=text size=8 name=similarityThreshold required value=0.5>"
        "<tr><td>Maximum connectivity<td class=centered><input type=text size=8 name=maximumConnectivity required value=3>"
        "<tr><td>Gene graph<td class=centered><input type=text size=8 name=geneGraphName required>"
        "</table>"
        "<input type=submit value='Create a new gene graph as described above'>"
        "</form>";
}



void ExpressionMatrix::exploreGeneGraph(
    const vector<string>& request,
    ostream& html, const BrowserInformation& browserInformation)
{
    // Get parameters from the request.
    string geneGraphName;
    getParameterValue(request, "geneGraphName", geneGraphName);
    if(geneGraphName.empty()) {
        html << "Gene graph name is missing.";
        return;
    }
    GeneGraph& geneGraph = getGeneGraph(geneGraphName);
    geneGraph.computeLayout();

    string coloringOption = "black";
    getParameterValue(request, "coloringOption", coloringOption);

    string metaDataName;
    getParameterValue(request, "metaDataName", metaDataName);


    GeneGraph::SvgParameters svgParameters;
    string showEdges;
    getParameterValue(request, "showEdges", showEdges);
    svgParameters.showEdges = (showEdges == "on");
    string showVertexLabels;
    getParameterValue(request, "showVertexLabels", showVertexLabels);
    svgParameters.showVertexLabels = (showVertexLabels == "on") && browserInformation.isChrome;
    getParameterValue(request, "svgSizePixels", svgParameters.svgSizePixels);
    getParameterValue(request, "xShift", svgParameters.xShift);
    getParameterValue(request, "yShift", svgParameters.yShift);
    getParameterValue(request, "zoomFactor", svgParameters.zoomFactor);
    getParameterValue(request, "vertexSizeFactor", svgParameters.vertexSizeFactor);
    getParameterValue(request, "edgeThicknessFactor", svgParameters.edgeThicknessFactor);

    string metaDataMeaning = "category";
    getParameterValue(request, "metaDataMeaning", metaDataMeaning);


    // Color the gene graph as requested.
    if(coloringOption == "byMetaData") {
        colorGeneGraphByMetaData(geneGraph, metaDataName);
    } else {
        colorGeneGraphBlack(geneGraph);
    }



    html << "<h2>Gene graph " << geneGraphName << "</h2>";





    // Write out the JavaScript and html to allow manipulating the svg.
    html << R"###(
<br>
To pan  the graphics, use the mouse to drag an empty location in the graph.
<br>
<form id=mouseWheelForm>
Mouse wheel controls: 
<input id=mouseWheelFormDefault type="radio" onclick='updateMouseWheelFunction(this)' name="mouseWheelFunction" value="zoom" checked=checked>zoom 
<input type="radio" onclick='updateMouseWheelFunction(this)' name="mouseWheelFunction" value="graphicsSize">graphics size 
<input type="radio" onclick='updateMouseWheelFunction(this)' name="mouseWheelFunction" value="vertexSize">vertex size 
<input type="radio" onclick='updateMouseWheelFunction(this)' name="mouseWheelFunction" value="edgeThickness">edge thickness
</form>



<script>
document.getElementById("mouseWheelFormDefault").checked = true;
// Function called when the user clicks one of the radio buttons
// that control mouse wheel function.
var mouseWheelFunction = "zoom";
function updateMouseWheelFunction(element) 
{
    if(element.checked) {
        mouseWheelFunction = element.value;
    }
}

// Turn off scrolling with the mouse wheel.
// Otherwise it interacts with this page's use of the mouse wheel.
// This is browser dependent. As coded, works 
// in Firefox 58.0 and Chrome 64.0.
function stopWheel(e){
    if(!e) { 
        e = window.event; 
    } 
    e.stopPropagation();
    e.preventDefault();
    e.cancelBubble = false;
    return false;
}
document.addEventListener('DOMMouseScroll', stopWheel, false);
document.addEventListener('wheel', stopWheel, false);

var mouseDown = false;
var xMouse = 0;
var yMouse = 0;


// Set the viewBox of the svg object based on the current values of the global variables.
function setViewBox()
{
    var svg = document.getElementById("svgObject");
    var xMin = xCenter-halfViewBoxSize;
    var yMin = yCenter-halfViewBoxSize;
    var viewBox = xMin + " " + yMin + " " + 2.*halfViewBoxSize + " " + 2.*halfViewBoxSize;
    svg.setAttribute("viewBox",  viewBox);
}


// Extract the delta in pixels from a wheel event.
// The constant is set in such a way that Firefox and Chrome
// behave in approximately the same way.
function extractWheelDelta(e) {
    if (e.deltaMode == WheelEvent.DOM_DELTA_LINE) {
        return e.deltaY * 27;
    } else {
        return e.deltaY;
    }
}

function mouseDownHandler(event) {
    // alert("mouseDownHandler");
    xMouse = event.clientX;
    yMouse = event.clientY;
    mouseDown = true;
}
function mouseMoveHandler(event) {
    if(mouseDown) {
        // alert("mouseMoveHandler");
        var xDelta = event.clientX - xMouse;
        var yDelta = event.clientY - yMouse;
        xCenter -= pixelSize*xDelta;
        yCenter -= pixelSize*yDelta;
        xShift += pixelSize*xDelta;
        yShift += pixelSize*yDelta;
        xMouse = event.clientX;
        yMouse = event.clientY;
        setViewBox();    
    }
}
function mouseUpHandler() {
    // alert("mouseUpHandler");
    mouseDown = false;
}



// Function called when the user moves the mouse wheel.
function handleMouseWheelEvent(e) 
{
    if(mouseWheelFunction == "zoom") {
        zoomHandler(e);
    } else if(mouseWheelFunction == "graphicsSize") {
        handleGraphicsResizeEvent(e);
    } else if(mouseWheelFunction == "vertexSize") {
        handleVertexResizeEvent(e);
    } else if(mouseWheelFunction == "edgeThickness") {
        handleEdgeThicknessEvent(e);
    }
}


// Function called when the user moves the mouse wheel to resize the vertices.
function handleVertexResizeEvent(e) {
    var delta = extractWheelDelta(e);
    var factor = Math.exp(-.001*delta);
    vertexSizeFactor *= factor;
    
    // Resize all the vertices.
    var vertices = document.getElementById("vertices").childNodes;
    for(i=0; i<vertices.length; i++) {
        vertex = vertices[i];
        if(vertex.tagName == "circle") {
            vertex.transform.baseVal.getItem(1).setScale(vertexSizeFactor, vertexSizeFactor); 
        }
    }
    
    // Resize all the vertex labels.
    var vertexLabels = document.getElementById("vertexLabels").childNodes;
    for(i=0; i<vertexLabels.length; i++) {
        label = vertexLabels[i];
        if(label.tagName == "text") {
            label.transform.baseVal.getItem(1).setScale(vertexSizeFactor, vertexSizeFactor); 
        }
    }
}



// Function called when the user moves the mouse wheel to change edge thickness.
function handleEdgeThicknessEvent(e) {
    var delta = extractWheelDelta(e);
    var factor = Math.exp(-.001*delta);
    edgeThicknessFactor *= factor;
    // alert(edgeThicknessFactor);
    
    
    // Change the thickness of all the edges.
    var edges = document.getElementById("edges").childNodes;
    for(i=0; i<edges.length; i++) {
        if(edges[i].tagName == "line") {
            edges[i].style.strokeWidth *= factor;
        }
    }
}



// Function called when the user moves the mouse wheel to zoom in or out.
function zoomHandler(e) {
    var delta = extractWheelDelta(e);
    var factor = Math.exp(.001*delta);
    zoomFactor /= factor;
    halfViewBoxSize *= factor;
    pixelSize *= factor;
    setViewBox();
}



// Function called when the user moves the mouse wheel to resize the svg graphics.
function handleGraphicsResizeEvent(e) {
    var delta = extractWheelDelta(e);
    var factor = Math.exp(-.001*delta);
    svgSizePixels *= factor;
    pixelSize /= factor;
    
    svg = document.getElementById("svgObject");
    svg.setAttribute("width", factor*svg.getAttribute("width"));
    svg.setAttribute("height", factor*svg.getAttribute("height"));
}


// Function called when the user clicks "Redraw graph."
function prepareColoringFormForSubmit()
{
    document.getElementById("svgSizePixels").value = svgSizePixels;
    document.getElementById("xShift").value = xShift;
    document.getElementById("yShift").value = yShift;
    document.getElementById("zoomFactor").value = zoomFactor;
    document.getElementById("vertexSizeFactor").value = vertexSizeFactor;
    document.getElementById("edgeThicknessFactor").value = edgeThicknessFactor;
}


function highlightGene()
{
    var geneName = document.getElementById("geneToHighlight").value;
    var circle = document.getElementById(geneName);
    if(!circle) {
        alert(geneName + " not found in gene graph");
    } else {
        circle.setAttribute("fill", "red");
        circle.setAttribute("stroke", "green");
    }
    
}
</script>
)###";


    // Div to contain the svg graphics and the form on the right.
    html << "<div>";


    // Write the svg.
    html <<
        "<div style='float:left;' "
        "onmousedown='mouseDownHandler(event);' "
        "onmouseup='mouseUpHandler(event);' "
        "onmousemove='mouseMoveHandler(event);' "
        "onwheel='handleMouseWheelEvent(event);'>";
    geneGraph.writeSvg(html, svgParameters, *this);
    html << "</div>";

    // Svg display parameters get written to the html in Javascript code
    html <<
        "<script>"
        "var svgSizePixels = " << svgParameters.svgSizePixels << ";"
        "var xShift = " << svgParameters.xShift << ";"
        "var yShift = " << svgParameters.yShift << ";"
        "var zoomFactor = " << svgParameters.zoomFactor << ";"
        "var vertexSizeFactor = " << svgParameters.vertexSizeFactor << ";"
        "var edgeThicknessFactor = " << svgParameters.edgeThicknessFactor << ";"
        "var xCenter = " << svgParameters.xCenter << ";"
        "var yCenter = " << svgParameters.yCenter << ";"
        "var halfViewBoxSize = " << svgParameters.halfViewBoxSize << ";"
        "var pixelSize = " << svgParameters.pixelSize << ";"
        "</script>";



    // Form on the right of the graph.
    html << "<div>";



    // Form for coloring the graph.
    html <<
        "<div>"
        "<form id=coloringForm onsubmit='prepareColoringFormForSubmit()'>"
        "<br><input type=radio name=coloringOption value=black";
    if(coloringOption == "black") {
        html << " checked=checked";
    }
    html <<
        ">Color black"
        "<br><input type=radio name=coloringOption value=byMetaData";
    if(coloringOption == "byMetaData") {
        html << " checked=checked";
    }
    html <<
        ">Color by meta data field ";
    set<string> selectedMetaData;
    if(!metaDataName.empty()) {
        selectedMetaData.insert(metaDataName);
    }
    writeGeneMetaDataSelection(html, "metaDataName", selectedMetaData, false);
    html <<
        "<br><input type=checkbox name=showEdges title='The graph displays faster if the edges are not displayed'";
    if(svgParameters.showEdges) {
        html << " checked=checked";
    }
    html << ">Display edges";
    if(browserInformation.isChrome) {
        html <<
            "<br><input type=checkbox name=showVertexLabels";
        if(svgParameters.showVertexLabels) {
            html << " checked=checked";
        }
        html << ">Display gene labels";
    } else {
        html << "<br>Gene labels are only supported with the Chrome browser.";
    }
    html <<
        "<input type=text hidden id=geneGraphName name=geneGraphName value='" << geneGraphName << "'>"
        "<input type=text hidden id=svgSizePixels name=svgSizePixels>"
        "<input type=text hidden id=xShift name=xShift>"
        "<input type=text hidden id=yShift name=yShift>"
        "<input type=text hidden id=zoomFactor name=zoomFactor>"
        "<input type=text hidden id=vertexSizeFactor name=vertexSizeFactor>"
        "<input type=text hidden id=edgeThicknessFactor name=edgeThicknessFactor>"
        "<p><button type=submit>Redraw graph</button>"
        " If browser behaves erratically, try clicking on \"Redraw graph\"."
        "</form>";


    // Form to highlight a gene.
    html <<
        "<button onclick='highlightGene()'>Highlight this gene:</button>"
        "<input id=geneToHighlight type=text>";



    // End of div containing the forms on the right.
    html << "</div>";

    // End of div containing the svg graphics and the form on the right.
    html << "</div>";

    if(svgParameters.showVertexLabels) {
        html <<
            "<p>If gene labels don't display correctly, "
            "a slight panning will often fix the problem. "
            "This is caused by browser bugs displaying text in svg graphics."
            "To pan, use the mouse to drag an empty location in the graph.";

    }
}



void ExpressionMatrix::createGeneGraph(const vector<string>& request, ostream& html)
{
    string geneSetName;
    getParameterValue(request, "geneSetName", geneSetName);
    if(geneSetName.empty()) {
        html << "Gene set name is missing.";
        return;
    }

    string similarGenePairsName;
    getParameterValue(request, "similarGenePairsName", similarGenePairsName);
    if(similarGenePairsName.empty()) {
        html << "Similar gene pairs name is missing.";
        return;
    }

    int maximumConnectivity = 3;
    getParameterValue(request, "maximumConnectivity", maximumConnectivity);
    double similarityThreshold = 0.5;
    getParameterValue(request, "similarityThreshold", similarityThreshold);

    string geneGraphName;
    getParameterValue(request, "geneGraphName", geneGraphName);
    if(geneGraphName.empty()) {
        html << "Gene graph name is missing.";
        return;
    }



    html << "<h1>Create gene graph " << geneGraphName << "</h1>";
    html << "<pre>";
    createGeneGraph(html, geneGraphName, geneSetName, similarGenePairsName,
        maximumConnectivity, similarityThreshold);
    html << "</pre>";
    html << "<p>Gene graph " << geneGraphName << " was created."
        "<p><form action=exploreGeneGraph>"
        "<input type=text hidden name=geneGraphName value=" << geneGraphName <<
        "><input type=submit value=Continue></form>";

}



void ExpressionMatrix::removeGeneGraph(const vector<string>& request, ostream& html)
{
    string geneGraphName;
    getParameterValue(request, "geneGraphName", geneGraphName);
    if(geneGraphName.empty()) {
        html << "Gene graph name is missing.";
        return;
    }
    removeGeneGraph(geneGraphName);

    html << "<p>Gene graph " << geneGraphName << " was removed."
        "<p><form action=exploreGeneGraphs>"
        "<input type=submit value=Continue></form>";

}



ostream& ExpressionMatrix::writeSimilarGenePairsSelection(
    ostream& html,
    const string& selectName) const
{
    vector<string> availableSimilarGenePairs;
    getAvailableSimilarGenePairs(availableSimilarGenePairs);

    html << "<select";
    html << " title='Select one'";
    html << " name=" << selectName << " style='vertical-align:text-top;'>";
    html << "<option value=''></option>";
    for(const string& name: availableSimilarGenePairs) {
        html << "<option value='" << name;
        html << "'>" << name << "</option>";
    }
    html << "</select>";

    return html;

}



void ExpressionMatrix::removeSimilarGenePairs(const vector<string>& request, ostream& html)
{
    string similarGenePairsName;
    getParameterValue(request, "similarGenePairsName", similarGenePairsName);

    removeSimilarGenePairs(similarGenePairsName);

    html <<
        "<p>Removed similar gene pairs " << similarGenePairsName << "."
        "<p><form action=similarGenePairs><input type=submit value=Continue></form>";
}
