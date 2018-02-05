#include "ExpressionMatrix.hpp"
#include "Lsh.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "SignatureGraph.hpp"
#include "color.hpp"

#include <boost/graph/iteration_macros.hpp>

#include "fstream.hpp"



// Check that a signature graph does not exist,
// and throw an exception if it does.
void ExpressionMatrix::checkSignatureGraphDoesNotExist(
    const string& signatureGraphName) const
{
    if(signatureGraphs.find(signatureGraphName) != signatureGraphs.end()) {
        throw runtime_error("Signature graph " + signatureGraphName + " already exists.");
    }

}

// Return a reference to a signature graph with a given name,
// and throw and exception if not found.
SignatureGraph& ExpressionMatrix::getSignatureGraph(const string& signatureGraphName)
{
    const auto it = signatureGraphs.find(signatureGraphName);
    if(it == signatureGraphs.end()) {
        throw runtime_error("Signature graph " + signatureGraphName + " does not exists.");
    }
    return *(it->second);
}


// All cells with the same signature are aggregated
// into a single vertex of the signature graph.
void ExpressionMatrix::createSignatureGraph(
    const string& signatureGraphName,
    const string& cellSetName,
    const string& lshName,
    size_t minCellCount)
{
    checkSignatureGraphDoesNotExist(signatureGraphName);

    // Locate the cell set and verify that it is not empty.
    const auto& it = cellSets.cellSets.find(cellSetName);
    if(it == cellSets.cellSets.end()) {
        throw runtime_error("Cell set " + cellSetName + " does not exist.");
    }
    const MemoryMapped::Vector<CellId>& cellSet = *(it->second);
    const CellId cellCount = CellId(cellSet.size());
    if(cellCount == 0) {
        throw runtime_error("Cell set " + cellSetName + " is empty.");
    }

    // Access the Lsh object that will do the computation.
    Lsh lsh(directoryName + "/Lsh-" + lshName);
    if(lsh.cellCount() != cellCount) {
        throw runtime_error("LSH object " + lshName +
            " has a number of cells inconsistent with cell set " + cellSetName + ".");
    }
    const size_t lshBitCount = lsh.lshCount();
    cout << "Number of LSH signature bits is " << lshBitCount << "." << endl;

    // Gather cells with the same signature.
    // These are cell ids local to the cell set we are using to create the signature graph.
    map<BitSetPointer, vector<CellId> > signatureMap;
    for(CellId cellId=0; cellId<cellCount; cellId++) {
        signatureMap[lsh.getSignature(cellId)].push_back(cellId);
    }
    cout << "Found " << signatureMap.size();
    cout << " populated signatures out of " << (1ULL<<lshBitCount);
    cout << " total possible signatures." << endl;

#if 0
    // Create a histogram containing how many signatures have a given number of cells.
    {
        cout << timestamp << "Creating signature histogram." << endl;
        vector<size_t> histogram;
        for(const auto& p: signatureMap) {
            const size_t size = p.second.size();
            if(size >= histogram.size()) {
                histogram.resize(size+1, 0);
            }
            ++(histogram[size]);
        }
        ofstream csvOut("SignatureHistogram.csv");
        csvOut << "Number of cells,Number of signatures,"
            "Total number of cells,Cumulative total number of cells\n";
        size_t sum = 0;
        for(size_t i=0; i<histogram.size(); i++) {
            const size_t frequency = histogram[i];
            if(frequency) {
                sum += frequency*i;
                csvOut << i << "," << frequency << "," << frequency*i  << "," << sum << "\n";
            }
        }
    }
#endif

    // Create the signature graph.
    const std::shared_ptr<SignatureGraph> signatureGraphPointer =
        std::make_shared<SignatureGraph>();
    signatureGraphs.insert(make_pair(signatureGraphName, signatureGraphPointer));
    SignatureGraph& signatureGraph = *signatureGraphPointer;

    // Create the vertices of the signature graph.
    cout << timestamp << "Creating vertices of the signature graph." << endl;
    typedef SignatureGraph::vertex_descriptor vertex_descriptor;
    for(const auto& p: signatureMap) {
        const BitSetPointer signature = p.first;
        const vector<CellId> cells = p.second;
        if(cells.size() < minCellCount) {
            continue;
        }
        const vertex_descriptor v = add_vertex(signatureGraph);
        signatureGraph.vertexMap.insert(make_pair(signature, v));
        SignatureGraphVertex& vertex = signatureGraph[v];
        vertex.signature = signature;
        vertex.localCellIds = cells;
        vertex.globalCellIds.reserve(cells.size());
        for(const CellId localCellId: cells) {
            vertex.globalCellIds.push_back(cellSet[localCellId]);
        }
    }
    cout << "The signature graph has " << num_vertices(signatureGraph) << " vertices." << endl;

    // Create the edges of the signature graph.
    cout << timestamp << "Creating edges of the signature graph." << endl;
    signatureGraph.createEdges(lshBitCount);
    cout << "The signature graph has " << num_edges(signatureGraph) << " edges." << endl;
    cout << "Average connectivity is " <<
        (2.*double(num_edges(signatureGraph)))/double(num_vertices(signatureGraph)) << endl;

    // Write out the signature graph in Graphviz format.
    // signatureGraph.writeGraphviz("SignatureGraph.dot");

    // Write out the signature graph in svg format.
    // This gives us more flexibility than using svg to create svg output.
    SignatureGraph::SvgParameters svgParameters;
    // svgParameters.hideEdges = true;
    signatureGraph.writeSvg("SignatureGraph.svg", svgParameters);

    cout << timestamp << "createSignatureGraph ends." << endl;
}


void ExpressionMatrix::removeSignatureGraph(const string& signatureGraphName)
{
    const auto it = signatureGraphs.find(signatureGraphName);
    if(it == signatureGraphs.end()) {
        throw runtime_error("Signature graph " + signatureGraphName + " does not exists.");
    }
    signatureGraphs.erase(it);
}


void ExpressionMatrix::colorBlack(SignatureGraph& graph) const
{
    const pair<string, double> p("black", 1.);
    BGL_FORALL_VERTICES(v, graph, SignatureGraph) {
        SignatureGraphVertex& vertex = graph[v];
        vertex.colors.resize(1);
        vertex.colors.front() = p;
    }

}



void ExpressionMatrix::colorRandom(SignatureGraph& graph) const
{
    // Objects used for random number generation.
    const int seed = 231;
    std::mt19937 randomGenerator(seed);
    std::uniform_real_distribution<double> distribution(0, 1.);

    // Loop over all vertices.
    string vertexColor;
    BGL_FORALL_VERTICES(v, graph, SignatureGraph) {
        SignatureGraphVertex& vertex = graph[v];
        const double red = distribution(randomGenerator);
        const double green = distribution(randomGenerator);
        const double blue = distribution(randomGenerator);
        vertexColor = color(red, green, blue);
        vertex.colors.resize(1);
        vertex.colors.front() = make_pair(vertexColor, 1.);
    }
}



void ExpressionMatrix::colorByMetaDataInterpretedAsCategory(SignatureGraph& graph, const string& metaDataName) const
{

    // Count the number of cells with each of the meta data values.
    // The map is keyed by meta data values.
    map<string, CellId> metaDataValuesTable;
    const StringId metaDataNameId = cellMetaDataNames(metaDataName);
    BGL_FORALL_VERTICES(v, graph, SignatureGraph) {
        const SignatureGraphVertex& vertex = graph[v];
        for(const CellId globalCellId: vertex.globalCellIds) {
            const string metaDataValue = getCellMetaData(globalCellId, metaDataNameId);
            const auto it = metaDataValuesTable.find(metaDataValue);
            if(it == metaDataValuesTable.end()) {
                metaDataValuesTable.insert(make_pair(metaDataValue, 1));
            } else {
                ++(it->second);
            }
        }
    }

    // Sort the meta data values by decreasing number of cells.
    vector< pair<string, CellId> > metaDataValues(metaDataValuesTable.begin(), metaDataValuesTable.end());
    sort(metaDataValues.begin(), metaDataValues.end(), OrderPairsBySecondGreater< pair<string, CellId> >());

    // Assign colors to meta data values.
    map<string, string> colorTable;
    for(size_t colorIndex=0; colorIndex<metaDataValues.size(); colorIndex++) {
        const string color = (colorIndex < 12) ? colorPalette1(colorIndex) : "black";
        colorTable.insert(make_pair(metaDataValues[colorIndex].first, color));
    }


    // Now we can color the vertices.
    BGL_FORALL_VERTICES(v, graph, SignatureGraph) {
        SignatureGraphVertex& vertex = graph[v];
        vertex.colors.clear();
        for(const CellId globalCellId: vertex.globalCellIds) {
            const string metaDataValue = getCellMetaData(globalCellId, metaDataNameId);
            const string color = colorTable[metaDataValue];
            bool done = false;
            for(auto& p: vertex.colors) {
                if(p.first == color) {
                    ++(p.second);
                    done = true;
                    break;
                }
            }
            if(!done) {
                vertex.colors.push_back(make_pair(color, 1.));
            }
        }

        // Sort with the most frequent ones first.
        sort(vertex.colors.begin(), vertex.colors.end(), OrderPairsBySecondGreater< pair<string, double> >());

        // Normalize so the sum is 1.
        double sum = 0;
        for(const auto& p: vertex.colors) {
            sum += p.second;
        }
        const double factor = 1./sum;
        for(auto& p: vertex.colors) {
            p.second *= factor;
        }
    }

}



void ExpressionMatrix::colorByMetaDataInterpretedAsColor(SignatureGraph&, const string& metaDataName) const
{
    throw runtime_error("Signature graph coloring by meta data interpreted as color is not implemented.");
}



void ExpressionMatrix::colorByMetaDataInterpretedAsNumber(SignatureGraph&, const string& metaDataName) const
{
    throw runtime_error("Signature graph coloring by meta data interpreted as number is not implemented.");
}

