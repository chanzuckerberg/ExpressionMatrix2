// The gene graph is a graph in which each vertex
// corresponds to a gene.
// An undirected edge is created between two vertices
// if the there is good similarity between the
// expression vectors of the corresponding genes.

#ifndef CZI_EXPRESSION_MATRIX2_GENE_GRAPH_HPP
#define CZI_EXPRESSION_MATRIX2_GENE_GRAPH_HPP

#include "Ids.hpp"

#include "boost_array.hpp"
#include <boost/graph/adjacency_list.hpp>

#include "iosfwd.hpp"
#include "string.hpp"



namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        class GeneGraph;
        class GeneGraphVertex;
        class GeneGraphEdge;
        class GeneSet;

        // The base class for class CellGraph.
        using GeneGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::undirectedS,
            GeneGraphVertex,
            GeneGraphEdge>;

    }
}



// A vertex of the gene graph.
class ChanZuckerberg::ExpressionMatrix2::GeneGraphVertex {
public:
    GeneId globalGeneId = invalidGeneId;    // Gene id as defined by the ExpressionMatrix.
    GeneId localGeneId = invalidGeneId;     // Gene id as defined by the GeneSet.
    array<double, 2> position;
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
        const GeneSet&,
        const string& similarGenePairsName,
        double similarityThreshold,
        size_t maxConnectivity
        );

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


    // Write out the gene graph in SVG format.
    void writeSvg(
        const string& fileName,
        SvgParameters&);
    void writeSvg(
        ostream& s,
        SvgParameters&);

};


#endif
