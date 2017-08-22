#ifndef CZI_EXPRESSION_MATRIX2_CLUSTER_GRAPH_HPP
#define CZI_EXPRESSION_MATRIX2_CLUSTER_GRAPH_HPP

// In the cluster graph, each vertex represents a cluster of the cell similarity graph.

#include "Ids.hpp"
#include <boost/graph/adjacency_list.hpp>

#include "iosfwd.hpp"
#include "string.hpp"
#include "vector.hpp"


namespace ChanZuckerberg {
	namespace ExpressionMatrix2 {

        class ClusterGraph;
        class ClusterGraphVertex;
        class ClusterGraphEdge;

        // The base class for class ClusterGraph.
        using ClusterGraphBaseClass = boost::adjacency_list<
            boost::setS,	// No parallel edges
            boost::listS,
            boost::undirectedS,
			ClusterGraphVertex,
			ClusterGraphEdge
			>;

        class CellSimilarityGraph;
        class GeneSet;

        namespace MemoryMapped {
        	template<class StringId> class StringTable;
        }
	}
}



class ChanZuckerberg::ExpressionMatrix2::ClusterGraphVertex {
public:

	// The clusterId of the cluster represented by this vertex.
	// This is the same as the clusterId for all the CellSimilarityGraph
	// vertices that correspond to this ClusterGraphVertex.
	uint32_t clusterId;

	// The cells in the cluster represented by this vertex.
	vector<CellId> cells;

	// The average gene expression for these cells.
	// This is a vector of size equal to the number of genes
	// in the gene set used to create the cell similarity graph.
	vector<double> averageGeneExpression;
};



class ChanZuckerberg::ExpressionMatrix2::ClusterGraphEdge {
public:
	double similarity;
};



class ChanZuckerberg::ExpressionMatrix2::ClusterGraph : public ClusterGraphBaseClass {
public:

	// Create the ClusterGraph from the CellSimilarityGraph.
	// This uses the clusterId stored in each CellSimilarityGraphVertex.
	ClusterGraph(const CellSimilarityGraph&);

	// Store in each edge the similarity of the two clusters, computed using the clusters
	// average expression stored in each vertex.
	void computeSimilarities();

	// Remove the vertices that correspond to small clusters.
	void removeSmallVertices(size_t clusterSizeThreshold);

	// Remove edges with low similarity.
	void removeWeakEdges(double similarityThreshold);

	// Make it a k-nn graph.
	// For each vertex, keep the best k edges.
	void makeKnn(size_t k);

	// Write in Graphviz format.
    void write(ostream&, const GeneSet&, const MemoryMapped::StringTable<GeneId>& geneNames) const;
    void write(const string& fileName, const GeneSet&, const MemoryMapped::StringTable<GeneId>& geneNames) const;

private:

    class Writer {
    public:
        Writer(const ClusterGraph&, const GeneSet&, const MemoryMapped::StringTable<GeneId>& geneNames);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const ClusterGraph& graph;
        const GeneSet& geneSet;
        const MemoryMapped::StringTable<GeneId>& geneNames;
    private:
        // Compute font size for a vertex  given number of cells.
        static int fontSize(size_t cellCount);
        // Compute font size for an edge  given numbers of cells of the two vertices.
        static int fontSize(size_t cellCount0, size_t cellCount1);
    };
};



#endif

