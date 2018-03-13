#include "ExpressionMatrix.hpp"
#include "SimilarPairs.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include <boost/algorithm/string.hpp>

#include <regex>



void ExpressionMatrix::removeGeneSet(
    const string& geneSetName)
{
    // Sanity check: prevent removal of the AllGenes gene set.
    if(geneSetName == "AllGenes") {
        throw runtime_error("Gene set AllGenes cannot be removed.");
    }

    // Locate the gene set to be removed.
    const auto it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    GeneSet& geneSet = it->second;

    // Remove it from disk.
    geneSet.remove();

    // Remove it from our table of gene sets.
    geneSets.erase(it);
}



// Return the gene set with a given name.
// Throw an exception if a gene set with this name does not exist.
GeneSet& ExpressionMatrix::getGeneSet(const string& geneSetName)
{
    const auto it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    return it->second;
}
const GeneSet& ExpressionMatrix::getGeneSet(const string& geneSetName) const
{
    const auto it = geneSets.find(geneSetName);
    if(it == geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " does not exist.");
    }
    return it->second;
}


// Get the global gene ids of the genes in a gene set.
vector<GeneId> ExpressionMatrix::getGeneSetGenes(const string& geneSetName)
{
    GeneSet& geneSet = getGeneSet(geneSetName);
    vector<GeneId> genes(geneSet.size());
    copy(geneSet.genes().begin(), geneSet.genes().end(), genes.begin());
    return genes;
}



// Create a new gene set consisting of genes whose name matches a given regular expression.
bool ExpressionMatrix::createGeneSetFromRegex(const string& geneSetName, const string& regexString)
{
    // Check if a gene set with this name already exists.
    if(geneSets.find(geneSetName) != geneSets.end()) {
        return false;
    }

    // Create the regular expression we are going to match.
    const std::regex regex(regexString);

    // Create the new gene set.
    GeneSet& geneSet = geneSets[geneSetName];
    geneSet.createNew(directoryName + "/GeneSet-" + geneSetName);
    for(GeneId geneId = 0; geneId != geneCount(); geneId++) {
        const string geneName = geneNames[geneId];
        if(std::regex_match(geneName, regex)) {
            geneSet.addGene(geneId);
        }
    }
    geneSet.sort();
    return true;
}


// Create a gene set consisting of the genes with names passed in a vector.
// Names that don't correspond to valid gene names are ignored.
// Returns true if successful, false if the specified gene set already exists.
bool ExpressionMatrix::createGeneSetFromGeneNames(
    const string& geneSetName,
    const vector<string>& geneNamesVector,
    int& ignoredCount,
    int& emptyCount)
{
    // Check if a gene set with this name already exists.
    if(geneSets.find(geneSetName) != geneSets.end()) {
        return false;
    }

    // Create the new gene set.
    GeneSet& geneSet = geneSets[geneSetName];
    geneSet.createNew(directoryName + "/GeneSet-" + geneSetName);
    ignoredCount = 0;
    emptyCount = 0;
    for(const string& geneName: geneNamesVector) {
        if(geneName.empty()) {
            ++emptyCount;
            continue;
        }
        const StringId stringId = geneNames(geneName);
        if(stringId == geneNames.invalidStringId) {
            ++ignoredCount;
        } else {
            geneSet.addGene(GeneId(stringId));
        }
    }
    geneSet.sort();

    return true;
}



// Same as above, but throw an exception if any of the gene names are empty
// or do not correspond to any gene.
void ExpressionMatrix::createGeneSetFromGeneNames(
    const string& geneSetName,
    const vector<string>& geneNamesVector)
{
    // Check if a gene set with this name already exists.
    if(geneSets.find(geneSetName) != geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " already exists.");
    }

    // Create the new gene set.
    GeneSet& geneSet = geneSets[geneSetName];
    geneSet.createNew(directoryName + "/GeneSet-" + geneSetName);
    for(const string& geneName: geneNamesVector) {
        if(geneName.empty()) {
            throw runtime_error("Empty gene name while creating gene set " + geneSetName);
        }
        const StringId stringId = geneNames(geneName);
        if(stringId == geneNames.invalidStringId) {
            throw runtime_error("Unknown gene name name " + geneName + " while creating gene set " + geneSetName);
        } else {
            geneSet.addGene(GeneId(stringId));
        }
    }
    geneSet.sort();

}



// Create a gene set consisting of the genes with ids passed in a vector.
void ExpressionMatrix::createGeneSetFromGeneIds(
    const string& geneSetName,
    const vector<GeneId>& geneIds)
{
    // Check if a gene set with this name already exists.
    if(geneSets.find(geneSetName) != geneSets.end()) {
        throw runtime_error("Gene set " + geneSetName + " already exists.");
    }

    // Create the new gene set.
    GeneSet& geneSet = geneSets[geneSetName];
    geneSet.createNew(directoryName + "/GeneSet-" + geneSetName);
    for(const GeneId& geneId: geneIds) {
        geneSet.addGene(geneId);
    }
    geneSet.sort();

}



bool ExpressionMatrix::createGeneSetIntersection(const string& inputSetsNames, const string& outputSetName)
{
    return createGeneSetIntersectionOrUnion(inputSetsNames, outputSetName, false);
}
bool ExpressionMatrix::createGeneSetUnion(const string& inputSetsNames, const string& outputSetName)
{
    return createGeneSetIntersectionOrUnion(inputSetsNames, outputSetName, true);
}
bool ExpressionMatrix::createGeneSetIntersectionOrUnion(
    const string& commaSeparatedInputSetsNames,
    const string& outputSetName,
    bool doUnion)
{
    // See if a gene set with the name of the output gene set already exists.
    if(geneSets.find(outputSetName) != geneSets.end()) {
        cout << "Gene set " << outputSetName << " already exists." << endl;
        return false;
    }

    // Parse the input gene sets.
    vector<string> inputSetsNames;
    boost::algorithm::split(inputSetsNames, commaSeparatedInputSetsNames, boost::is_any_of(","));

    // Check that all input gene sets exist.
    for(const string& inputSetName: inputSetsNames) {
        if(geneSets.find(inputSetName) == geneSets.end()) {
            cout << "gene set " << inputSetName << " does not exists." << endl;
            return false;
        }
    }

    // Compute the intersection or union.
    vector<GeneId> outputSetGenes;
    for(size_t i=0; i<inputSetsNames.size(); i++) {
        const string& inputSetName = inputSetsNames[i];
        vector<GeneId> inputSetGenes;
        geneSets[inputSetName].assertIsSorted();
        geneSets[inputSetName].getSortedGenes(inputSetGenes);
        if(i == 0) {
            outputSetGenes = inputSetGenes;
        } else {
            vector<GeneId> newOutputSetGenes;
            if(doUnion) {
                std::set_union(
                    outputSetGenes.begin(), outputSetGenes.end(),
                    inputSetGenes.begin(), inputSetGenes.end(),
                    back_inserter(newOutputSetGenes));
            } else {
                std::set_intersection(
                    outputSetGenes.begin(), outputSetGenes.end(),
                    inputSetGenes.begin(), inputSetGenes.end(),
                    back_inserter(newOutputSetGenes));
            }
            outputSetGenes.swap(newOutputSetGenes);
        }
    }



    // Store this gene set.
    GeneSet& outputGeneSet = geneSets[outputSetName];
    outputGeneSet.createNew(directoryName + "/GeneSet-" + outputSetName);
    for(const GeneId geneId: outputSetGenes) {
        outputGeneSet.addGene(geneId);
    }
    outputGeneSet.sort();
    return true;
}



bool ExpressionMatrix::createGeneSetDifference(
    const string& inputSetName0,
    const string& inputSetName1,
    const string& outputSetName)
{
    // See if a gene set with the name of the output gene set already exists.
    if(geneSets.find(outputSetName) != geneSets.end()) {
        cout << "Gene set " << outputSetName << " already exists." << endl;
        return false;
    }



    // Locate the input gene sets.
    const auto it0 = geneSets.find(inputSetName0);
    if(it0 == geneSets.end()) {
        cout << "Gene set " << inputSetName0 << " does not exists." << endl;
        return false;
    }
    GeneSet& inputSet0 = it0->second;
    inputSet0.assertIsSorted();
    const auto it1 = geneSets.find(inputSetName1);
    if(it1 == geneSets.end()) {
        cout << "Gene set " << inputSetName1 << " does not exists." << endl;
        return false;
    }
    GeneSet& inputSet1 = it1->second;
    inputSet1.assertIsSorted();


    // Compute the difference.
    vector<GeneId> inputSet0Genes;
    vector<GeneId> inputSet1Genes;
    inputSet0.getSortedGenes(inputSet0Genes);
    inputSet1.getSortedGenes(inputSet1Genes);
    vector<GeneId> outputSetGenes;
    std::set_difference(
        inputSet0Genes.begin(), inputSet0Genes.end(),
        inputSet1Genes.begin(), inputSet1Genes.end(),
        back_inserter(outputSetGenes));



    // Store this gene set.
    GeneSet& outputGeneSet = geneSets[outputSetName];
    outputGeneSet.createNew(directoryName + "/GeneSet-" + outputSetName);
    for(const GeneId geneId: outputSetGenes) {
        outputGeneSet.addGene(geneId);
    }
    outputGeneSet.sort();
    return true;
}


// Create a gene set consisting of genes that are expressed in at
// least a minimum number of cells.
// The resulting gene set will contain all of the
// genes in the specified existing gene set that are
// expressed in at least the specified minimum number of cells,
// counting only cells in the specified cell set.
void ExpressionMatrix::createWellExpressedGeneSet(
    const string& inputGeneSetName,
    const string& inputCellSetName,
    const string& outputGeneSetName,
    CellId minCellCount
    )
{
    // See if a gene set with the name of the output gene set already exists.
    if(geneSets.find(outputGeneSetName) != geneSets.end()) {
        throw runtime_error("Gene set " + outputGeneSetName + " already exists.");
    }

    // Access the input gene set and cell set.
    const GeneSet& inputGeneSet = getGeneSet(inputGeneSetName);
    inputGeneSet.assertIsSorted();
    const CellSet& inputCellSet = cellSet(inputCellSetName);

    // Counter to hold the number of cells that each gene
    // in our input gene set is expressed in.
    // Indexed by the local gene id in the input gene set.
    vector<GeneId> cellCount(inputGeneSet.size(), 0);

    // Loop over all cells in the given cell set.
    for(CellId localCellId=0; localCellId!=inputCellSet.size(); localCellId++) {
        const CellId globalCellId = inputCellSet[localCellId];
        const auto& expressionVector = cellExpressionCounts[globalCellId];

        // Loop over the expression vector of this cell.
        for(const auto& p: expressionVector) {
            const GeneId globalGeneId = p.first;
            const GeneId localGeneId = inputGeneSet.getLocalGeneId(globalGeneId);
            if(localGeneId != invalidGeneId) {
                CZI_ASSERT(localGeneId < inputGeneSet.size());
                ++(cellCount[localGeneId]);
            }
        }
    }

    // Create the output gene set.
    GeneSet& outputGeneSet = geneSets[outputGeneSetName];
    outputGeneSet.createNew(directoryName + "/GeneSet-" + outputGeneSetName);
    for(GeneId localGeneId=0; localGeneId<inputGeneSet.size(); localGeneId++) {
        if(cellCount[localGeneId] >= minCellCount) {
            const GeneId globalGeneId = inputGeneSet.getGlobalGeneId(localGeneId);
            outputGeneSet.addGene(globalGeneId);
        }
    }
    outputGeneSet.sort();
}



// Returns the names of the gene sets in the geneSets map that are identical
// to the gene set of a SimilarPairs object with the given name.
// Note that there could be zero, one, or multiple gene sets
// that satisfy this condition.
vector<string> ExpressionMatrix::geneSetNamesFromSimilarPairsName(const string& similarPairsName) const
{
    // Open the existing SimilarPairs object.
    const SimilarPairs similarPairs(directoryName, similarPairsName, true);

    // Start with no gene sets.
    vector<string> geneSetNames;

    // Loop over our map of gene sets.
    for(auto it=geneSets.begin(); it!=geneSets.end(); ++it) {
        if(it->second == similarPairs.getGeneSet()) {
            geneSetNames.push_back(it->first);
        }

    }

    // Return the names we found.
    return geneSetNames;
}
