#ifndef CZI_EXPRESSION_MATRIX2_SIGNATURE_GRAPH_HPP
#define CZI_EXPRESSION_MATRIX2_SIGNATURE_GRAPH_HPP

// In the signature graph, each vertex represents
// the set of all cells with a given signature.

#include "BitSet.hpp"
#include "Ids.hpp"

#include <boost/graph/adjacency_list.hpp>

#include "array.hpp"
#include "iosfwd.hpp"
#include "map.hpp"
#include "utility.hpp"



namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        class SignatureGraph;
        class SignatureGraphEdge;
        class SignatureGraphVertex;

        // The base class for class SignatureGraph.
        typedef boost::adjacency_list<
            boost::vecS,
            boost::vecS,
            boost::undirectedS,
            SignatureGraphVertex,
            SignatureGraphEdge> SignatureGraphBaseClass;
    }
}


class ChanZuckerberg::ExpressionMatrix2::SignatureGraphVertex {
public:

    // The signature common to all cells of this vertex.
    BitSetPointer signature;

    // The cells with this signature.
    vector<CellId> localCellIds;    // Local to cell set used to create the signature graph.
    vector<CellId> globalCellIds;
    CellId cellCount() const
    {
        return CellId(localCellIds.size());
    }

    // The position of this vertex in the 2-D graph layout.
    array<double, 2> position;

    // The colors to be used to render this vertex.
    // If more than onem, the vertex is rendered as a pie chart.
    // Each color comes with a weight that determines
    // its fraction of the pie chart.
    vector< pair<string, double> > colors;
};



class ChanZuckerberg::ExpressionMatrix2::SignatureGraphEdge {
public:
};



class ChanZuckerberg::ExpressionMatrix2::SignatureGraph :
    public SignatureGraphBaseClass {
public:

    // The vertex map gives the vertex with a given signature.
    map<BitSetPointer, vertex_descriptor> vertexMap;

    // Given the vertices, create the edges.
    void createEdges(size_t lshBitCount);

    // Write out the signature graph in Graphviz format.
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    // Parameters to control svg output.
    class SvgParameters {
    public:
        // Parameters that can be set arbitrarily.
        bool hideEdges = false;
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


    // Write out the signature graph in SVG format.
    void writeSvg(
        const string& fileName,
        SvgParameters&);
    void writeSvg(
        ostream& s,
        SvgParameters&);


private:
    class Writer {
    public:
        Writer(const SignatureGraph&);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const SignatureGraph& graph;
    };

    // Use Graphviz to compute the graph layout and store it in the vertex positions.
    void computeLayout();
    bool layoutWasComputed = false;
};



#endif
