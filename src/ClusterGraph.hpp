#ifndef CZI_EXPRESSION_MATRIX2_CLUSTER_GRAPH_HPP
#define CZI_EXPRESSION_MATRIX2_CLUSTER_GRAPH_HPP

// In the cluster graph, each vertex represents a cluster of the cell graph.

#include "Ids.hpp"
#include <boost/graph/adjacency_list.hpp>

#include "iosfwd.hpp"
#include "map.hpp"
#include "string.hpp"
#include "vector.hpp"


namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        class ClusterGraph;
        class ClusterGraphVertex;
        class ClusterGraphEdge;
        class ClusterGraphCreationParameters;

        // The base class for class ClusterGraph.
        using ClusterGraphBaseClass = boost::adjacency_list<
            boost::setS,    // No parallel edges
            boost::listS,
            boost::undirectedS,
            ClusterGraphVertex,
            ClusterGraphEdge
            >;

        class CellGraph;
        class GeneSet;

        namespace MemoryMapped {
            template<class StringId> class StringTable;
        }
    }
}



// Creation parameters for a ClusterGraph.
// These control the label propagation algorithm.
class ChanZuckerberg::ExpressionMatrix2::ClusterGraphCreationParameters {
public:
    size_t stableIterationCount = 3;    // Stop after this many iterations without changes.
    size_t maxIterationCount = 100;     // Stop after this many iterations no matter what.
    size_t seed = 231;                  // To initialize label propagation algorithm.
    size_t minClusterSize = 100;        // Minimum number of cells for a cluster to be retained.
    size_t maxConnectivity = 3;
    double similarityThreshold = 0.5;   // For edges of the cluster graph.
};



class ChanZuckerberg::ExpressionMatrix2::ClusterGraphVertex {
public:

	// The clusterId of the cluster represented by this vertex.
	// This is the same as the clusterId for all the CellGraph
	// vertices that correspond to this ClusterGraphVertex.
	uint32_t clusterId;

	// The cells in the cluster represented by this vertex.
	vector<CellId> cells;

	// The average gene expression for these cells.
	// This is a vector of size equal to the number of genes
	// in the gene set used to create the cell graph.
	vector<double> averageGeneExpression;
};



class ChanZuckerberg::ExpressionMatrix2::ClusterGraphEdge {
public:
	double similarity;
};



class ChanZuckerberg::ExpressionMatrix2::ClusterGraph : public ClusterGraphBaseClass {
public:

    // Create the ClusterGraph from the CellGraph.
    // This uses the clusterId stored in each CellGraphVertex.
    ClusterGraph(const CellGraph&, const GeneSet& geneSet);

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
    void write(
        ostream&,
        const string& clusterGraphName,
        const MemoryMapped::StringTable<GeneId>& geneNames) const;
    void write(
        const string& fileName,
        const string& clusterGraphName,
        const MemoryMapped::StringTable<GeneId>& geneNames) const;

    // Layout in svg and pdf format, stored in memory.
    string svg;
    string pdf;
    bool computeLayout(
        size_t timeoutSeconds,
        const string& clusterGraphName,
        const MemoryMapped::StringTable<GeneId>& geneNames);
    bool hasLayout() const
    {
        return (!svg.empty()) && (!pdf.empty());
    }

    // Maps clusterId to vertex_descriptor.
    map<uint32_t, vertex_descriptor> vertexMap;

private:

    vector<GeneId> geneSet;

    class Writer {
    public:
        Writer(
            const ClusterGraph&,
            const string& clusterGraphName,
            const vector<GeneId>&,
            const MemoryMapped::StringTable<GeneId>& geneNames);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const ClusterGraph& graph;
        const string& clusterGraphName;
        const vector<GeneId>& geneSet;
        const MemoryMapped::StringTable<GeneId>& geneNames;
    private:
        // Compute font size for a vertex  given number of cells.
        static int fontSize(size_t cellCount);
        // Compute font size for an edge  given numbers of cells of the two vertices.
        static int fontSize(size_t cellCount0, size_t cellCount1);
    };
};



#endif

