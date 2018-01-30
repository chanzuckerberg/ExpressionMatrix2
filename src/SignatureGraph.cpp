#include "SignatureGraph.hpp"
#include "color.hpp"
#include "filesystem.hpp"
#include "orderPairs.hpp"
using namespace ChanZuckerberg::ExpressionMatrix2;

#include "boost/algorithm/string.hpp"
#include <boost/graph/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>
#include "boost_lexical_cast.hpp"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "fstream.hpp"
#include "stdexcept.hpp"
#include "utility.hpp"



// Given the vertices, create the edges.
void SignatureGraph::createEdges(size_t lshBitCount)
{
    SignatureGraph& graph = *this;

    // Loop over all vertices.
    BitSet signature1(lshBitCount);
    BGL_FORALL_VERTICES(v0, graph, SignatureGraph) {
        const SignatureGraphVertex& vertex0 = graph[v0];
        const BitSetPointer signature0 = vertex0.signature;

        // For each zero bit, add an edge to the vertex
        // with the same bit set to 1.
        for(size_t bit=0; bit!=lshBitCount; bit++) {
            if(vertex0.signature.get(bit) == 0) {
                signature1 = signature0;
                signature1.set(bit);
                const auto it1 = vertexMap.find(signature1);
                if(it1 != vertexMap.end()) {
                    const vertex_descriptor v1 = it1->second;
                    add_edge(v0, v1, graph);
                }
            }
        }

    }
}



// Write out the signature graph in Graphviz format.
void SignatureGraph::writeGraphviz(const string& fileName) const
{
    ofstream s(fileName);
    writeGraphviz(s);
}
void SignatureGraph::writeGraphviz(ostream& s) const
{
    Writer writer(*this);
    boost::write_graphviz(s, *this, writer, writer, writer);
}
SignatureGraph::Writer::Writer(const SignatureGraph& graph) :
    graph(graph)
{
}



void SignatureGraph::Writer::operator()(std::ostream& s) const
{
    // s << "size=100 ratio=expand;\n";
    s << "node [shape=point];\n";
    // s << "edge [style=invis];n"; // len=10 penwidth=0.1
}


// Write out a vertex of the cell graph.
void SignatureGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const SignatureGraphVertex& vertex = graph[v];
    const double vertexSize = 1.e-2 * sqrt(double(vertex.cellCount));

    s << "[";
    s << "width=" << vertexSize;
    s << "]";
}



void SignatureGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = source(e, graph);
    const SignatureGraphVertex& vertex0 = graph[v0];
    const SignatureGraphVertex& vertex1 = graph[v0];
    double edgeWeight = vertex0.cellCount * vertex1.cellCount;

    s << "[weight=\"" << edgeWeight << "\"]";
}



// Use Graphviz to compute the graph layout and store it in the vertex positions.
void SignatureGraph::computeLayout()
{
    if(layoutWasComputed) {
        return;
    }

    using filesystem::remove;

    // Write the graph in Graphviz format.
    const string uuid = boost::uuids::to_string(boost::uuids::uuid(boost::uuids::random_generator()()));
    const string dotFileName = "/dev/shm/SignatureGraph-" + uuid + ".dot";
    writeGraphviz(dotFileName);

    // Run sfdp with output in Graphviz plain format.
    // See https://www.graphviz.org/doc/info/output.html#d:plain.
    const string dotPlainFileName = dotFileName + ".plain";
    const int systemReturnCode = ::system(("sfdp -O -Tplain " + dotFileName).c_str());
    const int sfdpReturnCode = WEXITSTATUS(systemReturnCode);   // Man page for system is not super clear on this.
    if(sfdpReturnCode!=0 && sfdpReturnCode!=1) {    // Sfdp returns 1 if built without triangulation library.
        remove(dotFileName);
        remove(dotPlainFileName);
        throw runtime_error("Error " +
            std::to_string(systemReturnCode) + " " +
            std::to_string(sfdpReturnCode) +
            " running sfdp.");
    }


    // Extract vertex positions from the output.
    ifstream file(dotPlainFileName);
    string line;
    vector<string> tokens;
    while(true) {

        // Get a line.
        getline(file, line);
        if(!file) {
            break;
        }

        // Parse it.
        boost::algorithm::split(tokens, line, boost::algorithm::is_any_of(" "));

        // Only parse lines that describe vertices.
        CZI_ASSERT(tokens.size() >= 1);
        if(tokens[0] != "node") {
            continue;
        }
        CZI_ASSERT(tokens.size() >= 4);

        // Extract the positions for this vertex.
        try {
            const vertex_descriptor v = lexical_cast<vertex_descriptor>(tokens[1]);
            SignatureGraphVertex& vertex = (*this)[v];
            vertex.position[0] = lexical_cast<double>(tokens[2]);
            vertex.position[1] = lexical_cast<double>(tokens[3]);
        } catch(std::exception& e) {
            const string message = "Error processing the following line of " + dotPlainFileName + ": " +  line + "\nError is: " + e.what();
            throw runtime_error(message);
        }
    }
    layoutWasComputed = true;

    // Remove the files we created.
    remove(dotFileName);
    remove(dotPlainFileName);
}



// Write out the signature graph in Graphviz format.
void SignatureGraph::writeSvg(
    const string& fileName,
    const SvgParameters& svgParameters)
{
    ofstream s(fileName);
    writeSvg(s, svgParameters);
}
void SignatureGraph::writeSvg(
    ostream& s,
    const SvgParameters& svgParameters)
{
    // Make sure the layout was computed.
    computeLayout();

    // Objects used for random number generation (to generate random colors for the vertices).
    const int seed = 231;
    std::mt19937 randomGenerator(seed);
    std::uniform_real_distribution<double> distribution(0, 1.);



    // Compute minimum and maximum of the vertices coordinates,
    // and the maximum number of cells.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = xMin;
    double yMax = xMax;
    CellId maxCellCount = 0;
    const SignatureGraph& graph = *this;
    BGL_FORALL_VERTICES(v, graph, SignatureGraph) {
        const SignatureGraphVertex& vertex = graph[v];
        const double x = vertex.position[0];
        const double y = vertex.position[1];
        xMin = min(xMin, x);
        xMax = max(xMax, x);
        yMin = min(yMin, y);
        yMax = max(yMax, y);
        maxCellCount = max(maxCellCount, vertex.cellCount);
    }

    // Center and size of the square bounding box containing all the vertices.
    const double xBoundingBoxCenter = (xMin + xMax) / 2.;
    const double yBoundingBoxCenter = (yMin + yMax) / 2.;
    const double xBoundingBoxSize = xMax - xMin;
    const double yBoundingBoxSize = yMax - yMin;
    const double boundingBoxSize = max(xBoundingBoxSize, yBoundingBoxSize);

    // Viewbox parameters.
    const double xViewBoxCenter = xBoundingBoxCenter - svgParameters.xShift;
    const double yViewBoxCenter = yBoundingBoxCenter - svgParameters.yShift;
    const double viewBoxSize = 1.05 * boundingBoxSize / svgParameters.zoomFactor;
    const double xMinViewBox = xViewBoxCenter - viewBoxSize/2.;
    const double yMinViewBox = yViewBoxCenter - viewBoxSize/2.;



    // Compute a table of vertices ordered by decreasing number of cells.
    vector< pair<vertex_descriptor, CellId> > sortedVertices;
    BGL_FORALL_VERTICES(v, graph, SignatureGraph) {
        const SignatureGraphVertex& vertex = graph[v];
        sortedVertices.push_back(make_pair(v, vertex.cellCount));
    }
    sort(sortedVertices.begin(), sortedVertices.end(), OrderPairsBySecondGreater< pair<vertex_descriptor, CellId> >());



    // Start the svg object.
    s <<
        "<svg "
        "width='" << svgParameters.svgSizePixels << "' "
        "height='" << svgParameters.svgSizePixels << "' "
        "viewBox='" << xMinViewBox << " " << yMinViewBox << " " << viewBoxSize << " " << viewBoxSize << "'"
        ">";

    // Draw the edges before the vertices, to avoid obscuring the vertices.
    if(!svgParameters.hideEdges) {
        const double unscaledEdgeThickness = 5.e-4 * boundingBoxSize;
        BGL_FORALL_EDGES(e, graph, SignatureGraph) {
            const vertex_descriptor v1 = source(e, graph);
            const vertex_descriptor v2 = target(e, graph);
            const SignatureGraphVertex& vertex1 = graph[v1];
            const SignatureGraphVertex& vertex2 = graph[v2];
            s << "<line x1='" << vertex1.position[0] << "' y1='" << vertex1.position[1] << "'";
            s << " x2='" << vertex2.position[0] << "' y2='" << vertex2.position[1] << "'";

            s << " style='stroke:black;stroke-width:" <<
                unscaledEdgeThickness * svgParameters.edgeThicknessFactor << "' />";
        }
    }



    // Write the vertices in order of decreasing size.
    // This way we mitigate obscuring of vertices by other vertices.
    const double largestVertexUnscaledRadius = 0.03 * boundingBoxSize;
    for(const auto& p: sortedVertices) {
        const vertex_descriptor v = p.first;
        const SignatureGraphVertex& vertex = graph[v];
        const double x = vertex.position[0];
        const double y = vertex.position[1];
        const double vertexRadius =
            svgParameters.vertexSizeFactor *
            largestVertexUnscaledRadius *
            sqrt(double(vertex.cellCount)/double(maxCellCount));
        const double red = distribution(randomGenerator);
        const double green = distribution(randomGenerator);
        const double blue = distribution(randomGenerator);
        const string vertexColor = color(red, green, blue);
        s << "<circle cx='" << x << "' cy='" << y << "' r='" <<
            vertexRadius << "' stroke='none' fill='" << vertexColor << "'></circle>";

    }



    // End the svg object.
    s << "</svg>";
}
