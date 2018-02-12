#include "GeneGraph.hpp"
#include "CZI_ASSERT.hpp"
using namespace ChanZuckerberg::ExpressionMatrix2;

#include "fstream.hpp"


GeneGraph::GeneGraph(
    const GeneSet&,
    const string& similarGenePairsName,
    double similarityThreshold,
    size_t maxConnectivity
    )
{
    CZI_ASSERT(0);
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
