// In the cluster graph, each vertex represents a cluster of the cell similarity graph.



#include "ClusterGraph.hpp"
#include "CellSimilarityGraph.hpp"
#include "CZI_ASSERT.hpp"
#include "deduplicate.hpp"
#include "MemoryMappedStringTable.hpp"
#include "orderPairs.hpp"
#include "regressionCoefficient.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/graphviz.hpp>

#include "algorithm.hpp"
#include "fstream.hpp"
#include "map.hpp"
#include "stdexcept.hpp"
#include "utility.hpp"



// Create the ClusterGraph from the CellSimilarityGraph.
// This uses the clusterId stored in each CellSimilarityGraphVertex.
ClusterGraph::ClusterGraph(const CellSimilarityGraph& cellSimilarityGraph)
{

	// Construct the vertices of the ClusterGraph.
	map<uint32_t, vertex_descriptor> vertexMap;	// Maps clusterId to vertex_descriptor.
	BGL_FORALL_VERTICES(cv, cellSimilarityGraph, CellSimilarityGraph) {
		const CellSimilarityGraphVertex& cVertex = cellSimilarityGraph[cv];
		const uint32_t clusterId = cVertex.clusterId;

		// Look for a vertex for this cluster.
		const auto it = vertexMap.find(clusterId);

		// If we don't have a vertex for this clusterId, create one now.
		if(it == vertexMap.end()) {
			vertex_descriptor v = add_vertex(*this);
			vertexMap.insert(make_pair(clusterId, v));
			ClusterGraphVertex& vertex = (*this)[v];
			vertex.clusterId = clusterId;
			vertex.cells.push_back(cVertex.cellId);
		}

		// If we already have a vertex for this clusterId, add this cell to that vertex.
		else {
			const vertex_descriptor v = it->second;
			ClusterGraphVertex& vertex = (*this)[v];
			CZI_ASSERT(vertex.clusterId == clusterId);
			vertex.cells.push_back(cVertex.cellId);
		}
	}


	// Create the edges by looping over all edges of the CellSimilarityGraph.
	BGL_FORALL_EDGES(ce, cellSimilarityGraph, CellSimilarityGraph) {

		// Find the CellSimilarityGraph vertices of this edge.
		const CellSimilarityGraph::vertex_descriptor cv0 = source(ce, cellSimilarityGraph);
		const CellSimilarityGraph::vertex_descriptor cv1 = target(ce, cellSimilarityGraph);
		const CellSimilarityGraphVertex& cVertex0 = cellSimilarityGraph[cv0];
		const CellSimilarityGraphVertex& cVertex1 = cellSimilarityGraph[cv1];

		// Find the corresponding vertices in the ClusterGraph.
		const auto it0 = vertexMap.find(cVertex0.clusterId);
		const auto it1 = vertexMap.find(cVertex1.clusterId);
		CZI_ASSERT(it0 != vertexMap.end());
		CZI_ASSERT(it1 != vertexMap.end());
		const vertex_descriptor v0 = it0->second;
		const vertex_descriptor v1 = it1->second;

		// If the vertices are distinct, add the edge. If the edge already exists, it will not be created,
		// because we use boost::setS for the edgeList template argument.
		if(v0 != v1) {
			add_edge(v0, v1, *this);
		}
	}


}



// Store in each edge the similarity of the two clusters, computed using the clusters
// average expression stored in each vertex.
void ClusterGraph::computeSimilarities()
{
	ClusterGraph& graph = *this;

	BGL_FORALL_EDGES(e, graph, ClusterGraph) {
		const vertex_descriptor v0 = source(e, graph);
		const vertex_descriptor v1 = target(e, graph);

		const ClusterGraphVertex& vertex0 = graph[v0];
		const ClusterGraphVertex& vertex1 = graph[v1];

		graph[e].similarity = regressionCoefficient(
			vertex0.averageGeneExpression,
			vertex1.averageGeneExpression);
	}
}


// Remove the vertices that correspond to small clusters.
void ClusterGraph::removeSmallVertices(size_t clusterSizeThreshold)
{
	vector<vertex_descriptor> verticesToBeRemoved;
	BGL_FORALL_VERTICES(cv, *this, ClusterGraph) {
		if((*this)[cv].cells.size() < clusterSizeThreshold) {
			verticesToBeRemoved.push_back(cv);
		}
	}
	for(const vertex_descriptor cv: verticesToBeRemoved) {
		clear_vertex(cv, *this);
		remove_vertex(cv, *this);
	}
}



// Remove edges with low similarity.
void ClusterGraph::removeWeakEdges(double similarityThreshold)
{
	vector<edge_descriptor> edgesToBeRemoved;
	BGL_FORALL_EDGES(e, *this, ClusterGraph) {
		if((*this)[e].similarity < similarityThreshold) {
			edgesToBeRemoved.push_back(e);
		}
	}

	for(const edge_descriptor e: edgesToBeRemoved) {
		boost::remove_edge(e, *this);
	}
}



// Turn the cluster graph into a k-nn graph.
// For each vertex, keep the best k edges.
void ClusterGraph::makeKnn(size_t k)
{
	// Find the edges to be kept.
	vector<edge_descriptor> edgesToBeKept;
	vector< pair<double, edge_descriptor> > vertexEdges;
	BGL_FORALL_VERTICES(v, *this, ClusterGraph) {

		// Gather all the edges of this vertex.
		vertexEdges.clear();
		BGL_FORALL_OUTEDGES(v, e, *this, ClusterGraph) {
			vertexEdges.push_back(make_pair((*this)[e].similarity, e));
		}

		// Sort them by decreasing similarity.
		sort(vertexEdges.begin(), vertexEdges.end(), std::greater< pair<double, edge_descriptor> >());

		// Keep up to k.
		if(vertexEdges.size() > k) {
			vertexEdges.resize(k);
		}

		// Mark them as edges to be kept.
		for(const auto& p: vertexEdges) {
			edgesToBeKept.push_back(p.second);
		}

	}
	deduplicate(edgesToBeKept);



	// Find the edges to be removed.
	vector<edge_descriptor> edgesToBeRemoved;
	BGL_FORALL_EDGES(e, *this, ClusterGraph) {
		if(!binary_search(edgesToBeKept.begin(), edgesToBeKept.end(), e)) {
			edgesToBeRemoved.push_back(e);
		}
	}

	// Now remove them.
	for(const edge_descriptor e: edgesToBeRemoved) {
		boost::remove_edge(e, *this);
	}

}



// Write the graph in Graphviz format.
void ClusterGraph::write(const string& fileName, const MemoryMapped::StringTable<GeneId>& geneNames) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, geneNames);
}
void ClusterGraph::write(ostream& s, const MemoryMapped::StringTable<GeneId>& geneNames) const
{
    Writer writer(*this, geneNames);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&ClusterGraphVertex::clusterId, *this));
}

ClusterGraph::Writer::Writer(const ClusterGraph& graph, const MemoryMapped::StringTable<GeneId>& geneNames) :
    graph(graph), geneNames(geneNames)
{
}



void ClusterGraph::Writer::operator()(std::ostream& s) const
{
    s << "tooltip=\"Cluster graph\";\n";
    s << "node [shape=circle];\n";
    s << "edge [fontsize=8];\n";
}


// Write out a vertex of the cluster graph.
// Some of the constants used herfe may need to be made configurable.
void ClusterGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    // Get the vertex.
    const ClusterGraphVertex& vertex = graph[v];

    // Get sorted expression counts.
    vector< pair<GeneId, double> > sortedExpressionCounts;
    for(GeneId geneId=0; geneId<vertex.averageGeneExpression.size(); geneId++) {
    	sortedExpressionCounts.push_back(make_pair(geneId, vertex.averageGeneExpression[geneId]));
    }
    sort(sortedExpressionCounts.begin(), sortedExpressionCounts.end(), OrderPairsBySecondGreater< pair<GeneId, double> >());

    // Begin vertex attributes.
    s << "[";



#if 0
    // HTML-style Label.
    s << "label=< <table border='0' cellpadding='0'>";
    s << "<tr><td align='left'><b>Cluster</b></td><td align='right'><b> " << vertex.clusterId << "</b></td></tr>";
    s << "<tr><td align='left'><b>Cells</b></td><td align='right'><b> " << vertex.cells.size() << "</b></td></tr>";
    const auto oldPrecision = s.precision(3);
    for(const auto& p: sortedExpressionCounts) {
    	if(p.second < 0.2) {
    		break;
    	}
    	s << "<tr><td align='left'>" << geneNames[p.first] << "</td><td align='right'> " << p.second << "</td></tr>";
    }
    s << "</table>";
    s << ">";
    s.precision(oldPrecision);
#endif

    // Label.
    s << "label=\"";
    s << "Cluster " << vertex.clusterId << "\\n";
    s << vertex.cells.size() << " cells\\n";
    const auto oldPrecision = s.precision(3);
    for(const auto& p: sortedExpressionCounts) {
    	if(p.second < 0.2) {
    		break;
    	}
    	s << geneNames[p.first] << " " << p.second << "\\n";
    }
    s << "\"";
    s.precision(oldPrecision);


    // Font size.
    s << " fontsize=" << fontSize(vertex.cells.size());

    // Vertex size.
    // s << " width=" << 0.2 * sqrt(double(vertex.cells.size()));

    // Tooltip.
    s << " tooltip=\"Cluster " << vertex.clusterId << "\"";

    // End vertex attributes.
    s << "]";
}



void ClusterGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    // Get the edge information we need.
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const ClusterGraphVertex& vertex0 = graph[v0];
    const ClusterGraphVertex& vertex1 = graph[v1];
    const ClusterGraphEdge& edge = graph[e];

    // Begin edge attributes.
    s << "[";

    const auto oldPrecision = s.precision(2);
    const auto oldOptions = s.setf(std::ios::fixed);
    s << "label=\"" << edge.similarity << "\"";
    s.precision(oldPrecision);
    s.setf(oldOptions);

    // Font size.
    s << " fontsize=" << fontSize(vertex0.cells.size(), vertex1.cells.size());

    // End edge attributes.
    s << "]";
}



// Compute font size for a vertex  given number of cells.
int ClusterGraph::Writer::fontSize(size_t cellCount)
{
	// Make it proportional to a power of the number of cells.
	int size = int(2.5 * pow(double(cellCount), 0.2));

	// Not too large.
	size = min(size, 24);

	// Not to small.
	size = max(size, 6);

	return size;
}



// Compute font size for an edge  given numbers of cells of the two vertices.
int ClusterGraph::Writer::fontSize(size_t cellCount0, size_t cellCount1)
{
	// Just return the smaller font size of the two vertices.
	return fontSize(min(cellCount0, cellCount1));
}


