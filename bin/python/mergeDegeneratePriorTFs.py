
networkFile = sys.argv[1]
outFileBase = sys.argv[2]

netIn = open(networkFile,'r')
tfTargDic = dict()
for line in itertools.islice(netIn,1,None): # skip line 1, as it's a header
	# networkFile -- three column tab-separated format: TF, target genes, interaction weight (w/ header)
	lineInf = line.strip('\n').split('\t')
	tfName = lineInf[0]
	target = lineInf[1]
	weight = lineInf[2]
	targetWeight = target + '*&$%||' + weight
	try:
		tfTargDic[tfName].append(targetWeight)
	except KeyError:
		tfTargDic[tfName] = [targetWeight,]
netIn.close()

# 2. Determine TF overlaps
tfNames = sorted(tfTargDic.keys())
totTfs = len(tfNames)
tfMergers = dict() # key = TF, item = list of TFs that should be merged
overlaps = dict() # key = TF1-TF2, item = overlap
tfTargNums = dict() # key = TF, item = number of targets
for tf1ind in range(0,totTfs):
	tf1 = tfNames[tf1ind]
	tf1targets = set(tfTargDic[tf1])
	numTf1Targs = len(tf1targets)
	tfTargNums[tf1] = numTf1Targs
	for tf2ind in range(tf1ind+1,totTfs):
		tf2 = tfNames[tf2ind]
		tf2targets = set(tfTargDic[tf2])
		numTf2Targs = len(tf2targets)	
		tfTargNums[tf2] = numTf2Targs
		overlap = tf2targets.intersection(tf1targets)
		overlapSize = len(overlap)
		overlaps[tf1 + '_' + tf2] = overlapSize
		overlaps[tf2 + '_' + tf1] = overlapSize
		# check to see if the targets/interaction signs are identical
		if numTf1Targs == numTf2Targs and numTf1Targs == overlapSize:
			if tf1 in tfMergers.keys():
				tfMergers[tf1].append(tf2)
			else:
				tfMergers[tf1] = [tf1,tf2]
			if tf2 in tfMergers.keys():
				tfMergers[tf2].append(tf1)
			else: 
				tfMergers[tf2] = [tf1,tf2]

# 3. Output
allMergedTfs = set(tfMergers.keys())
usedMergedTfs = list()	# keep tracks of used TFs, so we don't output mergers twice		
printedTfs = list()
overlapsToPrint = dict() # key = printTfName, item = overlaps with other TFs
netOutFile = outFileBase + '_merged_sp.tsv'
overlapsOutFile = outFileBase + '_overlaps.txt'
targetTotalsFile = outFileBase + '_targetTotals.txt'
mergedTfsFile = outFileBase + '_mergedTfs.txt'
netOut = open(netOutFile,'w')
netOut.write('Regulator\tTarget\tWeight\n')
overlapsOut = open(overlapsOutFile,'w')
targetTotalsOut = open(targetTotalsFile,'w')
mergedTfsOut = open(mergedTfsFile,'w')

for tf in tfNames:
	# if a TF has been merged and the compound TF has not already been printed,
	# add change its name to the compound TF name (and print it)
	printIt = 0
	if tf in allMergedTfs and tf not in usedMergedTfs:
		mergedTfs = tfMergers[tf]
		for mTf in mergedTfs:
			usedMergedTfs.append(mTf)
		if len(mergedTfs) > 2:	
			tfPrint = '_'.join(mergedTfs[0:2]) + '...'
		else:
			tfPrint = '_'.join(mergedTfs)
		mergedTfsOut.write(tfPrint + '\t' + ', '.join(mergedTfs) + '\n')
		printIt = 1
	elif tf not in allMergedTfs:
		tfPrint = tf
		printIt = 1
	if printIt:
		targetTotalsOut.write(tfPrint + '\t' + str(tfTargNums[tf]) + '\n')
		for targ in tfTargDic[tf]:
			netOut.write(tfPrint + '\t' + '\t'.join(targ.split('*&$%||')) + '\n')
		overlapsToPrint[tfPrint] = []
		for tf2 in printedTfs:
			# if either name has been merged, use only the first
			indTf1 = tfPrint.split('_')
			indTf2 = tf2.split('_')
			overlapsToPrint[tfPrint].append(str(overlaps[indTf1[0] + '_' + indTf2[0]]))
		printedTfs.append(tfPrint)
mergedTfsOut.close()
targetTotalsOut.close()
netOut.close()
print mergedTfsFile
print targetTotalsFile
print netOutFile
# print overlaps
overlapsOut.write('\t' + '\t'.join(printedTfs) + '\n')
for tf in printedTfs:
	printStuff = overlapsToPrint[tf]
	totSims = len(printStuff)
	zerosToPad = len(printedTfs) - totSims - 1
	if len(printStuff) > 0:
		overlapsOut.write(tf + '\t' + '\t'.join(printStuff) + '\t' + str(tfTargNums[tf.split('_')[0]]) + '\t0'*zerosToPad + '\n')
	else:
		overlapsOut.write(tf + '\t' + str(tfTargNums[tf.split('_')[0]]) + '\t0'*zerosToPad + '\n')
overlapsOut.close()
print overlapsOutFile		
