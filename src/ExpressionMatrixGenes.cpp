// This file contain the implementation of functionality
// of class ExpressionMatrix related to genes.

#include "ExpressionMatrix.hpp"
#include "CZI_ASSERT.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;


// Add a gene.
// Returns true if the gene was added, false if it was already present.
bool ExpressionMatrix::addGene(const string& geneName)
{
    CZI_ASSERT(geneSets.find("AllGenes") != geneSets.end());

    const StringId stringId = geneNames(geneName);
    if(stringId == geneNames.invalidStringId) {

        // This gene did not already exist. Add it to our container of gene names.
        const GeneId geneId = GeneId(geneNames[geneName]);

        // Make space for meta data for this gene.
        geneMetaData.push_back();

        // Set the GeneName meta data for this gene.
        setGeneMetaData(geneId, "GeneName", geneName);

        // Add this gene to the AllGenes gene set.
        geneSets["AllGenes"].addGene(geneId);
        geneSets["AllGenes"].forceSorted(); // We guarantee that it remains sorted.

        // Return true to indicate that the gene was added.
        return true;
    } else {

        // This gene already exists.
        return false;
    }
}



// Set a meta data (name, value) pair for a given gene.
// If the name already exists for that gene, the value is replaced.
void ExpressionMatrix::setGeneMetaData(const string& geneName, const string& name, const string& value)
{
    const GeneId geneId = geneIdFromName(geneName);
    if(geneId == invalidGeneId) {
        throw runtime_error("Gene " + geneName + " does not exist.");
    }
    setGeneMetaData(geneId, name, value);
}
void ExpressionMatrix::setGeneMetaData(GeneId geneId, const string& name, const string& value)
{
    const StringId nameId = geneMetaDataNames[name];
    const StringId valueId = geneMetaDataValues[value];
    setGeneMetaData(geneId, nameId, valueId);

}
void ExpressionMatrix::setGeneMetaData(GeneId geneId, StringId nameId, const string& value)
{
    const StringId valueId = geneMetaDataValues[value];
    setGeneMetaData(geneId, nameId, valueId);

}
void ExpressionMatrix::setGeneMetaData(GeneId geneId, StringId nameId, StringId valueId)
{
    // Scan the existing meta data for this gene, looking for this name.
    for(auto& p: geneMetaData[geneId]) {
        if(p.first == nameId) {
            p.second = valueId; // The name already exists. replace the value.
            return;
        }
    }

    // The name did not exist for this gene. Add this (name, value) pair.
    geneMetaData.push_back(geneId, make_pair(nameId, valueId));
    incrementGeneMetaDataNameUsageCount(nameId);

}

// Return the value of a specified meta data field for a given gene.
// Returns an empty string if the gene does not have the specified meta data field.
string ExpressionMatrix::getGeneMetaData(GeneId geneId, const string& name) const
{
    // Find the string id of the name.
    // If it does not exist, return an empty string.
    const StringId nameId = geneMetaDataNames(name);
    if(nameId == geneMetaDataNames.invalidStringId) {
        return "";
    }
    return getGeneMetaData(geneId, nameId);
}
string ExpressionMatrix::getGeneMetaData(GeneId geneId, StringId nameId) const
{

    // Scan the name/value pairs for this gene, looking for nameId.
    for(const auto& metaDataPair: geneMetaData[geneId]) {
        if(metaDataPair.first == nameId) {
            const StringId valueId = metaDataPair.second;
            if(valueId == geneMetaDataValues.invalidStringId) {
                return "";  // Should never happen, but just in case.
            } else {
                return geneMetaDataValues[valueId];
            }
        }
    }

    // We did not find it. Return an empty string.
    return "";
}



void ExpressionMatrix::incrementGeneMetaDataNameUsageCount(StringId nameId)
{
    if(geneMetaDataNamesUsageCount.size() <= nameId) {
        // This is a new name.
        CZI_ASSERT(geneMetaDataNamesUsageCount.size() == nameId);
        geneMetaDataNamesUsageCount.push_back(1);
    } else {

        // This is an existing name.
        ++(geneMetaDataNamesUsageCount[nameId]);
    }

}



void ExpressionMatrix::decrementGeneMetaDataNameUsageCount(StringId nameId)
{
    CZI_ASSERT(nameId < geneMetaDataNamesUsageCount.size());
    CZI_ASSERT(geneMetaDataNamesUsageCount[nameId] > 0);
    --(geneMetaDataNamesUsageCount[nameId]);
}


