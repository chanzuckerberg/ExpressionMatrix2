#!/usr/bin/python3
import argparse
import csv
import os
import shutil
# from ExpressionMatrix2 import ExpressionMatrix, NormalizationMethod

def get_data_size(expression, metadata, separator):
	"""
	Auto computes size of data based on expression and metadata files. 
	"""
	gene_count, cell_count = csv_dim(expression, separator)
	m_cell_count, m_category_count = csv_dim(metadata, separator)
	return {
		'geneCapacity': int(gene_count * 2.5), 
		'cellCapacity': int(cell_count * 4), 
		'cellMetaDataNameCapacity': int(m_cell_count * 4), 
		'cellMetaDataValueCapacity': int(m_category_count * m_cell_count * 4)
	}

def csv_dim(filename, separator):
	"""
	Computes rows and columns of a csv file
	"""
	row = 0
	col = 0
	with open(filename) as fi:
		reader = csv.reader(fi, delimiter=separator)
		col = len(next(reader))
		row = sum(1 for row in fi) + 1 # first line consumed by csv reader
	return row, col

def create_expression_matrix(expression, metadata, separator, sizes, output):
	"""
	Creates expression matrix object from delimited files
	"""
	if os.path.exists(output):
		shutil.rmtree(output)
	e = ExpressionMatrix(
		directoryName = output,
		geneCapacity = sizes['geneCapacity'],
		cellCapacity = sizes['cellCapacity'],
		cellMetaDataNameCapacity = sizes['cellMetaDataNameCapacity'],
		cellMetaDataValueCapacity = sizes['cellMetaDataValueCapacity']
	)
	e.addCells(
		expressionCountsFileName = expression,
		expressionCountsFileSeparators = separator,
		cellMetaDataFileName = metadata,
		cellMetaDataFileSeparators = separator
	)
	return e


if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Creates data files for Expression Matrix 2 from csv data. "
		"The expression file must have cells as columns and genes as rows. By default the separator is ',' but you can override with the --separator option. "
		"The metadata file must have metadata categories as columns and cells as row. Cell id must be the first column. "
		"The output will be stored in a folder called data (this can be overridden with the --output-folder options. ")
	parser.add_argument('-e', '--expression-file', dest='expression_file', action='store', required=True, help='CSV file with expression data')
	parser.add_argument('-m', '--metadata-file', dest='metadata_file', action='store', required=True, help='CSV file with metadata for cells')
	parser.add_argument('--separator', dest='separator', action='store', default=',', help='Separator for csv files')
	parser.add_argument('--output-folder', dest='output', action='store', default='data', help='Output folder for EM2 data')
	parser.add_argument('--similar-pairs-name', dest='similarpairsall', action='store', default='LshAll', help='Name for similar pairs for all genes')
	parser.add_argument('--similar-pairs-highinfo', dest='similarpairshighinfo', action='store', default='LshHighInformationGenes', help='Name for similar pairs for high info genes')
	parser.add_argument('--highinfo-genes-name', dest='highinfogenes', action='store', default='HighInformationGenes', help='Name for high information genes set')
	parser.add_argument('--cluster-name', dest='clustername', action='store', default='EM2Cluster', help='Metadata key for clusters')
	parser.add_argument('--skip-clustering', dest='skipcluster', action='store_true',  help="Don't cluster data (if you data is too large to cluster all of it)")

	args = parser.parse_args()

	sizes = get_data_size(args.expression_file, args.metadata_file, args.separator)
	# Create expression matrix with sizes
	e = create_expression_matrix(args.expression_file, args.metadata_file, args.separator, sizes, args.output)
	# find similar pairs
	# name- lsh-all
	e.findSimilarPairs4(similarPairsName=args.similarpairsall)
	# create high info genes
	e.createGeneSetUsingInformationContent(normalizationMethod=NormalizationMethod.L2, geneInformationContentThreshold=2, newGeneSetName=args.highinfogenes)
	# find similar pairs for high info
	e.findSimilarPairs4(similarPairsName=args.similarpairshighinfo, geneSetName=args.highinfogenes)
	# create graph (not layout)
	# cluster and store in metadata (name option)
	if not args.skipcluster:
		e.createCellGraph('cellgraph', similarPairsName=args.similarpairshighinfo)
		e.createClusterGraph(cellGraphName='cellgraph', clusterGraphName='clustergraph')
		e.createMetaDataFromClusterGraph(clusterGraphName='clustergraph', metaDataName=args.clustername)




