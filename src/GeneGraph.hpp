// The gene graph is a graph in which each vertex
// corresponds to a gene.
// An undirected edge is created between two vertices
// if the there is good similarity between the
// expression vectors of the corresponding genes.

#ifndef CZI_EXPRESSION_MATRIX2_GENE_GRAPH_HPP
#define CZI_EXPRESSION_MATRIX2_GENE_GRAPH_HPP

#include "Ids.hpp"

#include <boost/graph/adjacency_list.hpp>

#include "array.hpp"
#include "iosfwd.hpp"
#include "map.hpp"
#include "string.hpp"



namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        class GeneGraph;
        class GeneGraphVertex;
        class GeneGraphEdge;
        class ExpressionMatrix;
        class GeneSet;

        // The base class for class CellGraph.
        using GeneGraphBaseClass = boost::adjacency_list<
            boost::setS,    // Prevent parallel edges.
            boost::listS,
            boost::undirectedS,
            GeneGraphVertex,
            GeneGraphEdge>;

    }
}



// A vertex of the gene graph.
// Each vertex corresponds to a gene.
class ChanZuckerberg::ExpressionMatrix2::GeneGraphVertex {
public:

    // The GeneId local to the gene set that was
    // used to create the gene graph.
    // This is not necessarily the same as the gene set
    // used to create the SimilarGenePairs object
    // used for graph creaton.
    // GeneId localGeneId = invalidGeneId;

    // The global GeneId known to the ExpressionMatrix object.
    GeneId globalGeneId = invalidGeneId;

    // Position of this vertex (gene) in the graph layout.
    array<double, 2> position;

    // Color to be used to display this vertex.
    string color;

    // Constructors.
    GeneGraphVertex() {}
    GeneGraphVertex(GeneId globalGeneId) :
        globalGeneId(globalGeneId) {}
};



// An edge of the cell graph.
class ChanZuckerberg::ExpressionMatrix2::GeneGraphEdge {
public:
    float similarity = -1.;

    GeneGraphEdge()
    {
    }

    GeneGraphEdge(float similarity) :
        similarity(similarity)
    {
    }
};



class ChanZuckerberg::ExpressionMatrix2::GeneGraph : public GeneGraphBaseClass {
public:

    // Use the constructors of the base class.
    using GeneGraphBaseClass::GeneGraphBaseClass;

    using Graph = GeneGraph ;
    Graph& graph()
    {
        return *this;
    }
    const Graph& graph() const
    {
        return *this;
    }

    GeneGraph(
        ostream& out,
        GeneSet&,
        const string& similarGenePairsName,
        double similarityThreshold,
        size_t maxConnectivity
        );

    // Write out the signature graph in Graphviz format.
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    // Parameters to control svg output.
    class SvgParameters {
    public:
        // Parameters that can be set arbitrarily.
        bool showEdges = true;
        bool showVertexLabels = true;
        double svgSizePixels = 500;
        double xShift = 0.;
        double yShift = 0.;
        double zoomFactor = 1.;
        double vertexSizeFactor = 1.;
        double edgeThicknessFactor = 1.;
        // Parameters derived from the above values.
        // These are filled in by writeSvg.
        double xCenter;
        double yCenter;
        double halfViewBoxSize;
        double pixelSize;
    };

    // Use Graphviz to compute the graph layout and store it in the vertex positions.
    void computeLayout();

    // Write out the gene graph in SVG format.
    void writeSvg(
        const string& fileName,
        SvgParameters&,
        const ExpressionMatrix&);
    void writeSvg(
        ostream& s,
        SvgParameters&,
        const ExpressionMatrix&);

private:
    // Map from global gene id to vertex descriptor.
    map<GeneId, vertex_descriptor> vertexTable;

    class Writer {
    public:
        Writer(const GeneGraph&);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const GeneGraph& graph;
    };

    bool layoutWasComputed = false;
};


#endif
