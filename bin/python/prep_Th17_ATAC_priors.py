#!/usr/bin/python

"""
prep_Th17_ATAC_priors.py -- take Dayanne's priors looking at different p-value cutoffs for 
motifs and peak-gene association rules, prepare for TRN inference
1.  priorsTable2Sparse.py
	convert a prior from table to sparse output
2.  mergeDegeneratePriorTFs.py
	Given degeneracy of TF motifs, many priors based on ATAC-seq data contain TFs with identical
	target genes and interaction strengths.  This function creates meta-TFs (whose name will 
	be a compound of individual TFs in the merge (e.g., TF1_TF2...)).  Note: because mergers
	often contained many TF names, we cut off at the first two, and then user can track down
	the rest of the group in the merger file below  
	USAGE:
		python ${scriptHome}/mergeDegeneratePriorTFs.py networkFile outFileBase
	INPUTS:
		networkFile -- three column tab-separated format: TF, target genes, interaction weight (w/ header)
		outFileBase -- output file base name
	OUTPUTS:
		0. new networkFile, where meta-TF names replace individual TFs with identical targets
		netOutFile = outFileBase + '_merged_sp.tsv'
		1. a target overlap table (TF X TF) w/ # of of overlapping targets
		overlapsOutFile = outFileBase + '_overlaps.txt'
		2. a table, containing number of targets per TF	
		targetTotalsFile = outFileBase + '_targetTotals.txt'
		3. a table, column 1 = abbreviated merged TF name, column 2 = all TFs included in merger
		mergedTfsFile = outFileBased + '_mergedTfs.txt'
3.  priorsSparse2Table.py
	convert a prior in a sparse, 3-column format (regulator, row, interaction strength) to a 
	table format (columns = regulators, rows = gene targets)
	Usage: python ${scriptHome}/priorsSparse2Table.py priorSparseFile priorTableOutFile quantCut
	INPUTS:
		priorSparseFile -- sparse, 3-column output: 
			regulator, target, and interaction strength
		priorTableOutFile -- name for the tabular prior
		quantCut -- value s.t. interactions with |quantile rank| > quantCut will be included
			in resulting interaction file		
	OUTPUTS:
		priorTableOutFile -- tab-delimited file, columns = regulators, rows = gene targets, 
			values = interaction strengths

"""	

import sys
import os
from shutil import copyfile

scriptHome = '/Users/caz3so/workspaces/infTRN_lassoStARS/priorParsingFxns'
	
priorDirectory = '/Users/caz3so/workspaces/infTRN_lassoStARS/Th17_example/inputs/priors/Tareian_prior'

outDir = '/Users/caz3so/workspaces/infTRN_lassoStARS/Th17_example/outputs'

# file base names for priors
fileBases = ['Th2_48hr_prior_count_binary',
	'Th2_48hr_prior_count',
	'Th17_48hr_prior_count_binary',
	'Th17_48hr_prior_count']

for fb0 in fileBases:
	fb = outDir + '/' + fb0
	copyfile(priorDirectory + '/' + fb0 + '.tsv',fb + '.tsv') 
	# convert table prior to sparse format
	sys.argv = [scriptHome + '/priorsTable2Sparse.py',
		fb + '.tsv',
		fb + '_sp.tsv',
		0]
	execfile(scriptHome + '/priorsTable2Sparse.py')

	##  run merged degenerate TFs
	sys.argv = [scriptHome + '/mergeDegeneratePriorTFs.py',
		fb + '_sp.tsv',
		fb ]
	execfile(scriptHome + '/mergeDegeneratePriorTFs.py')
	print " "
	
## get a matrix form network
for fb0 in fileBases:
	fb = outDir + '/' + fb0
	print fb	
	sys.argv = [scriptHome + '/priorsSparse2Table.py',
		fb + '_merged_sp.tsv',
		fb + '_merged.tsv',
		0]
	execfile(scriptHome + '/priorsSparse2Table.py')
	print " "



	