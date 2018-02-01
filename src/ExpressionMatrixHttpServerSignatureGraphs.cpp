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
    // Get parameters from the request.
    string signatureGraphName;
    getParameterValue(request, "signatureGraphName", signatureGraphName);
    if(signatureGraphName.empty()) {
        html << "Signature graph name is missing.";
        return;
    }
    SignatureGraph& signatureGraph = getSignatureGraph(signatureGraphName);

    string coloringOption = "random";
    getParameterValue(request, "coloringOption", coloringOption);

    string metaDataName;
    getParameterValue(request, "metaDataName", metaDataName);

    SignatureGraph::SvgParameters svgParameters;
    string hideEdges;
    getParameterValue(request, "hideEdges", hideEdges);
    svgParameters.hideEdges = (hideEdges == "on");
    getParameterValue(request, "svgSizePixels", svgParameters.svgSizePixels);
    getParameterValue(request, "xShift", svgParameters.xShift);
    getParameterValue(request, "yShift", svgParameters.yShift);
    getParameterValue(request, "zoomFactor", svgParameters.zoomFactor);
    getParameterValue(request, "vertexSizeFactor", svgParameters.vertexSizeFactor);
    getParameterValue(request, "edgeThicknessFactor", svgParameters.edgeThicknessFactor);



    html << "<h2>Signature graph " << signatureGraphName << "</h2>";


    // Color the signature graph as requested.
    if(coloringOption == "byMetaData") {
        colorByMetaDataInterpretedAsCategory(signatureGraph, metaDataName);
    } else if(coloringOption == "random") {
        colorRandom(signatureGraph);
    } else {
        colorBlack(signatureGraph);
    }



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
            // vertex.setAttribute("r", factor*vertex.getAttribute("r"));
            vertex.transform.baseVal.getItem(1).setScale(vertexSizeFactor, vertexSizeFactor); 
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
</script>
)###";


    // Div to contain the svg graphics and the form for coloring.
    html << "<div>";


    // Write the svg.
    html <<
        "<div style='float:left;'"
        "onmousedown='mouseDownHandler(event);'"
        "onmouseup='mouseUpHandler(event);'"
        "onmousemove='mouseMoveHandler(event);'"
        "onwheel='handleMouseWheelEvent(event);'>";
    signatureGraph.writeSvg(html, svgParameters);
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
        "<br><input type=radio name=coloringOption value=random";
    if(coloringOption == "random") {
        html << " checked=checked";
    }
    html <<
        ">Color randomly"
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
    writeMetaDataSelection(html, "metaDataName", selectedMetaData, false);
    html << "<br><input type=checkbox name=hideEdges";
    if(svgParameters.hideEdges) {
        html << " checked=checked";
    }
    html <<
        ">Hide edges"
        "<input type=text hidden id=signatureGraphName name=signatureGraphName value='" << signatureGraphName << "'>"
        "<input type=text hidden id=svgSizePixels name=svgSizePixels>"
        "<input type=text hidden id=xShift name=xShift>"
        "<input type=text hidden id=yShift name=yShift>"
        "<input type=text hidden id=zoomFactor name=zoomFactor>"
        "<input type=text hidden id=vertexSizeFactor name=vertexSizeFactor>"
        "<input type=text hidden id=edgeThicknessFactor name=edgeThicknessFactor>"
        "<p><button type=submit>Redraw graph</button>"
        "</form>"
        "</div>";



    // End the div containing the svg graphics and the form for coloring.
    html << "</div>";
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
