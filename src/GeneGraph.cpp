// CZI.
#include "GeneGraph.hpp"
#include "CZI_ASSERT.hpp"
#include "GeneSet.hpp"
#include "SimilarGenePairs.hpp"
using namespace ChanZuckerberg::ExpressionMatrix2;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard libraries.
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
    for(const GeneId localGeneId: geneSet.genes()) {
        const GeneId globalGeneId = geneSet.getGlobalGeneId(localGeneId);
        const vertex_descriptor v = add_vertex(GeneGraphVertex(localGeneId, globalGeneId), graph);
        vertexTable.insert(make_pair(localGeneId, v));
    }



    // Create the edges.
    for(const GeneId localGeneId0: geneSet.genes()) {
        const vertex_descriptor v0 = vertexTable[localGeneId0];

        // Only add the first maxConnectivity neighbors that are also in the gene set.
        // The similar genes are stored by decreasing similarity.
        size_t connectivity = 0;
        for(const auto& p: similarGenePairs[localGeneId0]) {
            const float similarity = p.second;
            if(similarity < similarityThreshold) {
                break;  // They are sorted by decreasing similarity, so no need to look at the rest.
            }
            const GeneId localGeneId1 = p.first;
            const auto it1 = vertexTable.find(localGeneId1);
            if(it1 != vertexTable.end()) {
                const vertex_descriptor v1 = vertexTable[localGeneId1];
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
        vertexTable.erase(graph[v].localGeneId);
        clear_vertex(v, graph);
        remove_vertex(v, graph);
    }

    out << "The gene graph has " << boost::num_vertices(graph);
    out << " vertices and " << num_edges(graph) << " edges\nafter";
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
        boost::get(&GeneGraphVertex::localGeneId, *this));
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
    const double vertexSize = 0.1;

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
    SvgParameters& svgParameters)
{
    ofstream file(fileName);
    writeSvg(file, svgParameters);
}
void GeneGraph::writeSvg(
    ostream& s,
    SvgParameters& svgParameters)
{
    CZI_ASSERT(0);
}
