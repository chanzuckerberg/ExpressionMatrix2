#include "SignatureGraph.hpp"
using namespace ChanZuckerberg::ExpressionMatrix2;

#include <boost/graph/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "fstream.hpp"



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
void SignatureGraph::write(const string& fileName) const
{
    ofstream s(fileName);
    write(s);
}
void SignatureGraph::write(ostream& s) const
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
    // s << "node [shape=point];\n";
    // s << "edge [style=invis];n"; // len=10 penwidth=0.1
}


// Write out a vertex of the cell graph.
void SignatureGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const SignatureGraphVertex& vertex = graph[v];
    const double vertexSize = 1.e-3 * sqrt(double(vertex.cellCount));

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
    double edgeWeight = vertex0.cellCount + vertex1.cellCount;
    edgeWeight *= edgeWeight;

    s << "[weight=\"" << edgeWeight << "\"]";
}



