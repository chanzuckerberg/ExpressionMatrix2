// The cell similarity graph is a graph in which each vertex
// corresponds to a cell.
// An undirected edge is created between two vertices
// if the there is good similarity between the
// expression vectors of the corresponding cells.

#ifndef CZI_EXPRESSION_MATRIX2_CELL_SIMILARITY_GRAPH_H
#define CZI_EXPRESSION_MATRIX2_CELL_SIMILARITY_GRAPH_H

#include "Ids.hpp"

#include "boost_array.hpp"
#include <boost/graph/adjacency_list.hpp>

#include "iosfwd.hpp"
#include "map.hpp"
#include "string.hpp"
#include "vector.hpp"



namespace ChanZuckerberg {
	namespace ExpressionMatrix2 {

        class CellSimilarityGraph;
        class CellSimilarityGraphVertex;
        class CellSimilarityGraphEdge;

        // The base class for classa CellSimilarityGraph.
        typedef boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::undirectedS,
            CellSimilarityGraphVertex,
            CellSimilarityGraphEdge> CellSimilarityGraphBaseClass;

        namespace MemoryMapped {
            template<class T> class Vector;
        }
	}
}



// A vertex of the cell similarity graph.
class ChanZuckerberg::ExpressionMatrix2::CellSimilarityGraphVertex {
public:
	CellId cellId = invalidCellId;
	CellSimilarityGraphVertex() {}
	CellSimilarityGraphVertex(CellId cellId) :
	    cellId(cellId) {}

    array<double, 2> position;
    uint32_t group = 0;
    uint32_t clusterId = 0;
	string color;
    double value;
};



// An edge of the cell similarity graph.
class ChanZuckerberg::ExpressionMatrix2::CellSimilarityGraphEdge {
public:
	float similarity = -1.;

	CellSimilarityGraphEdge() {}

	CellSimilarityGraphEdge(float similarity) :
		similarity(similarity) {}

	string color;
};



class ChanZuckerberg::ExpressionMatrix2::CellSimilarityGraph : public CellSimilarityGraphBaseClass {
public:
	using CellSimilarityGraphBaseClass::CellSimilarityGraphBaseClass;
	typedef CellSimilarityGraph Graph;
	Graph& graph() {return *this;}
	const Graph& graph() const {return *this;}

	CellSimilarityGraph(
	    const MemoryMapped::Vector<CellId>& cellSet, // The cell set to be used.
        const string& similarPairsName,              // The name of the SimilarPairs object to be used to create the graph.
        double similarityThreshold,                  // The minimum similarity to create an edge.
        size_t maxConnectivity                       // The maximum number of neighbors (k of the k-NN graph).
    );

	// Only keep an edge if it is one of the best k edges for either
	// of the two vertices. This turns the graph into a k-nearest-neighbor graph.
	void keepBestEdgesOnly(std::size_t k);

	// Write in Graphviz format.
	void write(ostream&) const;
	void write(const string& fileName) const;

	// Simple graph statistics.
	ostream& writeStatistics(ostream&) const;

	// Remove isolated vertices and returns\ the number of vertices that were removed
	size_t removeIsolatedVertices();

	// Use Graphviz to compute the graph layout and store it in the vertex positions.
	void computeLayout();

	// Clustering using the label propagation algorithm.
	// The cluster each vertex is assigned to is stored in the clusterId data member of the vertex.
	void labelPropagationClustering(
	    ostream&,
	    size_t seed,                            // Seed for random number generator.
	    size_t stableIterationCountThreshold,   // Stop after this many iterations without changes.
	    size_t maxIterationCount                // Stop after this many iterations no matter what.
	    );

	// Vertex table, keyed by cell id.
	map<CellId, vertex_descriptor> vertexTable;

	// Compute minimum and maximum coordinates of all the vertices.
	void computeCoordinateRange(
        double& xMin,
        double& xMax,
        double& yMin,
        double& yMax) const;

	// Assign integer colors to groups.
	// The same color can be used for multiple groups, but if two
	// groups are joined by one or more edges they must have distinct colors.
	// On return, colorTable[group] contains the integer color assigned to each group.
	// This processes the groups in increasing order beginning at group 0,
	// so it is best if the group numbers are all contiguous, starting at zero,
	// and in decreasing size of group.
	void assignColorsToGroups(vector<uint32_t>& colorTable);

	// Write the graph in svg format.
	// This does not use Graphviz. It uses the graph layout stored in the vertices,
	// and previously computed using Graphviz.
	// The last argument specifies the color assigned to each vertex group.
	// If empty, vertex groups are not used, and each vertex is drawn
	// with its own color.
	void writeSvg(
	    ostream& s,
	    double svgSizePixels,
	    double xViewBoxCenter,
	    double yViewBoxCenter,
	    double viewBoxHalfSize,
	    double vertexRadius,
	    double edgeThickness,
	    const map<int, string>& groupColors
	    ) const;

	class Writer {
	public:
		Writer(const Graph&);
		void operator()(ostream&) const;
		void operator()(ostream&, vertex_descriptor) const;
		void operator()(ostream&, edge_descriptor) const;
		const Graph& graph;
        double minEdgeSimilarity;
	};
};

#endif


