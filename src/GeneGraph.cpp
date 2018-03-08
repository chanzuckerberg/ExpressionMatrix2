// CZI.
#include "GeneGraph.hpp"
#include "CZI_ASSERT.hpp"
#include "GeneSet.hpp"
#include "ExpressionMatrix.hpp"
#include "SimilarGenePairs.hpp"
using namespace ChanZuckerberg::ExpressionMatrix2;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard libraries.
#include "algorithm.hpp"
#include "fstream.hpp"
#include "utility.hpp"


GeneGraph::GeneGraph(
    ostream& out,
    GeneSet& geneSet,
    const string& similarGenePairsName,
    double similarityThreshold,
    size_t maxConnectivity
    )
{
    GeneGraph& graph = *this;

    // Access the similar gene pairs.
    SimilarGenePairs similarGenePairs(similarGenePairsName, true);

    // Create the vertices.
    for(const GeneId globalGeneId: geneSet.genes()) {
        const vertex_descriptor v = add_vertex(GeneGraphVertex(globalGeneId), graph);
        vertexTable.insert(make_pair(globalGeneId, v));
    }



    // Create the edges.
    // Note that there are two gene sets involved: the gene set to
    // be used for graph creation and the gene set that was used
    // to create the SimilarPairs. If we want to make sure not to
    // lose edges, the former must be a subset of the latter.
    // However, for flexibility we do not check for this.
    for(const GeneId globalGeneId0: geneSet.genes()) {
        const vertex_descriptor v0 = vertexTable[globalGeneId0];

        // Find the local gene id (in the cell set of the SimilarGenePairs object)
        // corresponding to this global gene id.
        // If the gene set of the SimilarPairs object does not contain this gene,
        // we skip this gene.
        // This could result in missing some edges in the graph.
        const CellId localGeneId0 = similarGenePairs.getGeneSet().getLocalGeneId(globalGeneId0);
        if(localGeneId0 == invalidGeneId) {
            continue;
        }

        // Only add the first maxConnectivity neighbors that are also in the gene set.
        // The similar genes are stored by decreasing similarity.
        size_t connectivity = 0;
        for(const auto& p: similarGenePairs[localGeneId0]) {
            const float similarity = p.second;
            if(similarity < similarityThreshold) {
                break;  // They are sorted by decreasing similarity, so no need to look at the rest.
            }
            const GeneId localGeneId1 = p.first;
            const GeneId globalGeneId1 = similarGenePairs.getGeneSet().getGlobalGeneId(localGeneId1);
            const auto it1 = vertexTable.find(globalGeneId1);
            if(it1 != vertexTable.end()) {
                const vertex_descriptor v1 = it1->second;
                add_edge(v0, v1, GeneGraphEdge(similarity), graph);
                ++connectivity;
                if(connectivity == maxConnectivity) {
                    break;
                }
            }
        }
    }

    // Remove isolated vertices.
    vector<vertex_descriptor> isolatedVertices;
    BGL_FORALL_VERTICES(v, graph, GeneGraph) {
        if(out_degree(v, graph) == 0) {
            isolatedVertices.push_back(v);
        }
    }
    for(const vertex_descriptor v: isolatedVertices) {
        vertexTable.erase(graph[v].globalGeneId);
        clear_vertex(v, graph);
        remove_vertex(v, graph);
    }

    out << "The gene graph has " << boost::num_vertices(graph);
    out << " vertices and " << num_edges(graph) << " edges\nafter ";
    out << isolatedVertices.size() << " vertices were removed. "<< endl;
}

// Write out the signature graph in Graphviz format.
void GeneGraph::writeGraphviz(const string& fileName) const
{
    ofstream s(fileName);
    writeGraphviz(s);
}
void GeneGraph::writeGraphviz(ostream& s) const
{
    Writer writer(*this);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&GeneGraphVertex::globalGeneId, *this));
}
GeneGraph::Writer::Writer(const GeneGraph& graph) :
    graph(graph)
{
}



void GeneGraph::Writer::operator()(std::ostream& s) const
{
    // s << "size=100 ratio=expand;\n";
    s << "node [shape=point];\n";
    // s << "edge [style=invis];n"; // len=10 penwidth=0.1
}


// Write out a vertex of the cell graph.
void GeneGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const double vertexSize = 0.01;

    s << "[";
    s << "width=" << vertexSize;
    s << "]";
}



void GeneGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
}



// Use Graphviz to compute the graph layout and store it in the vertex positions.
void GeneGraph::computeLayout()
{
    if(layoutWasComputed) {
        return;
    }

    using filesystem::remove;

    // Write the graph in Graphviz format.
    const string uuid = boost::uuids::to_string(boost::uuids::uuid(boost::uuids::random_generator()()));
    const string dotFileName = "/dev/shm/GeneGraph-" + uuid + ".dot";
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
            const vertex_descriptor v = vertexTable[lexical_cast<GeneId>(tokens[1])];
            GeneGraphVertex& vertex = (*this)[v];
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


// Write out the signature graph in SVG format.
void GeneGraph::writeSvg(
    const string& fileName,
    SvgParameters& svgParameters,
    const ExpressionMatrix& expressionMatrix)
{
    ofstream file(fileName);
    writeSvg(file, svgParameters, expressionMatrix);
}
void GeneGraph::writeSvg(
    ostream& s,
    SvgParameters& svgParameters,
    const ExpressionMatrix& expressionMatrix)
{
    // Make sure the layout was computed.
    computeLayout();



    // Compute minimum and maximum of the vertices coordinates.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = xMin;
    double yMax = xMax;
    const GeneGraph& graph = *this;
    BGL_FORALL_VERTICES(v, graph, GeneGraph) {
        const GeneGraphVertex& vertex = graph[v];
        const double x = vertex.position[0];
        const double y = vertex.position[1];
        xMin = min(xMin, x);
        xMax = max(xMax, x);
        yMin = min(yMin, y);
        yMax = max(yMax, y);
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

    // Fill in the rest of the SvgParameters.
    svgParameters.xCenter = xViewBoxCenter;
    svgParameters.yCenter = yViewBoxCenter;
    svgParameters.halfViewBoxSize = viewBoxSize / 2.;
    svgParameters.pixelSize = viewBoxSize / svgParameters.svgSizePixels;



    // Start the svg object.
    s <<
        "<svg id=svgObject "
        "style='border:solid DarkBlue;margin:2;' "
        "width='" << svgParameters.svgSizePixels << "' "
        "height='" << svgParameters.svgSizePixels << "' "
        "viewBox='" << xMinViewBox << " " << yMinViewBox << " " << viewBoxSize << " " << viewBoxSize << "'"
        ">";


    // Draw the edges before the vertices, to avoid obscuring the vertices.
    if(svgParameters.showEdges) {
        s << "<g id=edges>";
        const double unscaledEdgeThickness = 5.e-4 * boundingBoxSize;
        BGL_FORALL_EDGES(e, graph, GeneGraph) {
            const vertex_descriptor v1 = source(e, graph);
            const vertex_descriptor v2 = target(e, graph);
            const GeneGraphVertex& vertex1 = graph[v1];
            const GeneGraphVertex& vertex2 = graph[v2];
            s << "<line x1='" << vertex1.position[0] << "' y1='" << vertex1.position[1] << "'";
            s << " x2='" << vertex2.position[0] << "' y2='" << vertex2.position[1] << "'";

            s << " style='stroke:black;stroke-width:" <<
                unscaledEdgeThickness * svgParameters.edgeThicknessFactor << "' />";
        }
        s << "</g>";
    }



    // Write the vertices.
    s << "<g id=vertices>"; // fill-opacity='0.5'
    const double vertexUnscaledRadius = 2.e-3 * boundingBoxSize;
    BGL_FORALL_VERTICES(v, graph, GeneGraph) {
        const GeneGraphVertex& vertex = graph[v];
        const double x = vertex.position[0];
        const double y = vertex.position[1];
        const double vertexRadius = vertexUnscaledRadius;
        const string geneName = expressionMatrix.geneName(vertex.globalGeneId);


        // Draw the vertex as a circle.
        s << "<circle cx='0' cy='0' r='" <<
            vertexRadius << "' stroke='none' fill='" << vertex.color << "'"
            " transform='translate(" << x << " " << y << ") scale(" << svgParameters.vertexSizeFactor << ")'"
            // " onclick='window.location=\"gene?geneId=" << geneName << "\";'"
            " onclick='window.open(\"gene?geneId=" << geneName << "\");'"
            " cursor=pointer id='" << geneName << "'><title>" << geneName << "</title></circle>";


    }
    s << "</g>";



    // Write the vertex labels.
    if(svgParameters.showVertexLabels) {
        s << "<g id=vertexLabels>"; // fill-opacity='0.5'
        BGL_FORALL_VERTICES(v, graph, GeneGraph) {
            const GeneGraphVertex& vertex = graph[v];
            const double x = vertex.position[0];
            const double y = vertex.position[1];
            const double vertexRadius = vertexUnscaledRadius;
            const string geneName = expressionMatrix.geneName(vertex.globalGeneId);
            const double labelOffset = 0.7*vertexRadius;

            // Add the label.
            s << "<text x='" << labelOffset << "' y='-" << labelOffset << "'" <<
                " transform='translate(" << x << " " << y << ") scale(" << svgParameters.vertexSizeFactor << ")'"
                " font-size='" << 2.*vertexRadius << "' font-weight='bold' fill='green'>" <<
                geneName << "</text>";

        }
        s << "</g>";
    }

    // End the svg object.
    s << "</svg>";
}
