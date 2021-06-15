#!/bin/python

# By Hagen Tilgner 02.2021
# Translated from awk to python3 by Anoushka Joglekar 02.2021
""" Takes in a list of mapped reads,
10X barcode-UMI pairs, and an ONT fastq file.
Identifies 28 mers on forward and reverse strand
with atleast 1 matching 11 mer in a 22 bp window.
Outputs a file which is then filtered for the minimum
edit distance per read on fwd and reverse strand
"""
import os
import sys
import argparse
import pandas as pd
import numpy as np
import Levenshtein as lv
from traceback import print_exc
import gzip
import re
import multiprocessing as mp
import itertools

def revComp(my_seq):
	base_comp = {'A':'T', 'C':'G','G':'C', 'T':'A', 'N':'N', " ":" "}
	lms = list(my_seq[:-1]) ## parse string into list of components
	lms.reverse()
	try:
		lms = [base_comp[base] for base in lms]
	except TypeError:
		pass
	lms = ''.join(lms)
	return(lms)

def processAI(inFile,cols):
	""" Read in 10X barcode-umi-gene file
	and convert to a DF b.c. that's
	easier to work with """
	readList = [[re.split('\t|\.',a.strip('\n'))[i] for i in cols] for a in open(inFile,'r')]
	truncSeq = [a[1]+a[2][0:6] for a in readList]
	bc_umi = [a[1]+a[2] for a in readList]
	df = pd.DataFrame(readList, columns = ['Gene','Barcode','UMI'])
	df['truncSeq'] = truncSeq
	df['bc_umi'] = bc_umi
	return(df)

def findElements(localSeq, tc, localPosition, kLen, bcUMI_len):
	""" Given a search space, use a sliding window to get k-mers,
	find perfect matches for the k-mer in the truncated 10X whitelist,
	return 28-mers (BC+UMI) which start with the corresponding 11-mers,
	also return the indices of the truncated BC+UMI sequence so we can get
	the 28-mers after exiting the function"""
	startScan = localPosition - 50
	if startScan < 0:
		startScan = 0
	endScan = localPosition + 50
	if endScan > len(localSeq)-kLen:
		endScan = len(localSeq)-kLen
	str_list = [localSeq[i:i+kLen] for i in range(startScan,endScan)]	## List of kmers
	B = [[(trunc[i:i+kLen],trunc) for i in range(kLen +1)] for trunc in tc]
	C = dict(itertools.chain.from_iterable(B))
	try:
		matching = [tc.index(C.get(i)) for i in str_list if C.get(i)!=None]
		if matching:
			match_str, matching = secondLevelMatching(localSeq, matching, str_list, tc)
		else:
			matching = 'NONE'
			match_str = 'NONE'
		return matching, match_str
	except ValueError:
		matching = 'NONE'
		match_str = 'NONE'
		return matching, match_str

def secondLevelMatching(localSeq, matching, str_list, tc):
	matches = [tc[m] for m in matching]
	subList = dict([(m,[m[i:i+11] for i in range(12)]) for m in matches])
	s_int_m = [list(set(str_list).intersection(subList[key])) for key in subList.keys()]
	substr_ix = list([[s_int_m[k][0],subList[list(subList.keys())[k]].index(s_int_m[k][0])] for k in range(len(subList))])
	ix = [localSeq.index(substr_ix[k][0])-substr_ix[k][1] for k in range(len(substr_ix))]
	match_str = [localSeq[x:x+28] for x in ix if x>=0]
	toDrop = [x for x in range(len(ix)) if ix[x] < 0]
	matching = list(np.delete(matching, toDrop))
	if not matching:
		matching = 'NONE'
		match_str = 'NONE'
	return match_str, matching

def strandAgnostic(localSeq, tc, localPosition, subDF):
	""" Have to look on forward and reverse strand. Also need to
	use index returned by findElements function to recover the
	untruncated BC+UMI 28-mer"""
	m, ms = findElements_old(localSeq, tc, localPosition, 11, 28)
	if m != 'NONE':
		bcUMI = [list(subDF['bc_umi'])[i] for i in m]
		l = len(bcUMI)
	else:
		bcUMI = 'NONE'
		l = 0
	return bcUMI, ms, l

def getMatches(rn, line, tc, colNames, gene, d, subDF):
	"""get matches on both forward and reverse strand
	and add to dict"""
	fwdSeq = line[0:bcPosFwd+100]
	revSeq = revComp(line[-(bcPosRev+100):])
	bcUMI_f, msf, lf = strandAgnostic(fwdSeq, tc, bcPosFwd, subDF)
	bcUMI_r, msr, lr = strandAgnostic(revSeq, tc, bcPosRev, subDF)
	outLine = [rn, gene, lf, bcUMI_f, msf, lr, bcUMI_r, msr]
	d1 = [(colNames[i],outLine[i]) for i in range(len(outLine))]
	for key, val in d1:
		d[key].append(val)
	return(d)


def chunkAndProcess(args, txdf, read2gene, d, legalGenes):
	""" Create child processes and allot chunks of read from
	fastq file per child.
	Read fastq file line by line. If read is in the allInfo
	file and gene is in the 10X barcode-umi-gene whitelist, then
	look for 11-mer matches in the two ends (forward and reverse strand)
	of the sequence"""
	step = 4	## read lines in file with a step size of 4
	startLine = 0
	sim_processes = args.numProc
	print("Number of simultaneous processes: ", sim_processes)
	readCount = args.readNum // 4
	toFork = (readCount // sim_processes) + 1
        # Set childNo to 1 to give first child that childNo.
	childNo = 1
	while readCount >= 1:
		isChild = os.fork()
		if isChild == 0:
			break
		else:
			if childNo == sim_processes:
				break
			else:
				startLine = startLine + (toFork * 4)
				readCount = readCount - toFork
				childNo += 1
	if isChild != 0:
		os.waitpid(isChild, 0)
		sys.exit()
	with gzip.open(args.fq,"rt",encoding ='utf-8') as f:
		for _ in zip(range(startLine), f):
			pass
		for lineno, line in enumerate(f, start = startLine):
			if lineno < startLine + (toFork * 4):
				if lineno % step == 0:
					rn = line.split(' ')[0][1:]
					if rn in read2gene:
						gene = read2gene[rn]
					else:
						gene = ''
						continue
				if lineno % step == 1 and len(line)>=401 and gene in legalGenes:
					line = line.strip('\n')
					subDF = txdf[txdf['Gene']==gene]
					tc = list(subDF['truncSeq'])
					if tc:
						d = getMatches(rn, line, tc, colNames, gene, d, subDF)
			else:
				break
	outDF = pd.DataFrame(data=d)
	outDF.to_csv("%s_%d" %(args.outName, childNo), sep = "\t",index=False, header=False)
	return

def main():
	""" Parse args, create an empty dictionary. Read in the files
	and start analysis"""
	args = parse_args()
	global colNames, bcPosFwd, bcPosRev, legalGenes
	d = {'ReadName':[], 'Gene':[], 'numFwd' : [],'tx_fwd':[], 'ont_fwd':[], 'numRev' : [],'tx_rev':[],'ont_rev':[]}
	colNames = list(d.keys())
	bcPosFwd = args.fwd
	bcPosRev = args.rev
	readList = [a.decode('utf8').strip('\n').split('\t') for a in gzip.open(args.allInfo,'rb')]
	read2gene = {a[0]:a[1].split(".")[0] for a in readList}
	legalGenes = [gene.strip('\n').split(".")[0] for gene in open(args.geneList,'r')]
	del readList
	print("Reading in 10x whitelist")
	txdf = processAI(args.tx,[1,3,4])
	print("Done, moving on to finding barcodes")
	chunkAndProcess(args, txdf, read2gene, d, legalGenes)
	return


def parse_args():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--tx', type=str, help='Gene-BC-UMI whitelist from 10x data')
	parser.add_argument('--allInfo', type=str, help='AllInfo file containing Read-Gene mapping information')
	parser.add_argument('--fq', type=str, help='ONT fastq.gz file')
	parser.add_argument('--geneList', type=str, help='List of allowed genes')
	parser.add_argument('--fwd', type=int, default=50, help='Expected position of BC on fwd strand')
	parser.add_argument('--rev', type=int, default=35, help='Expected position of BC on rev strand')
	parser.add_argument('--outName', type=str, help='name of output file')
	parser.add_argument('--numProc', type=int, default=12, help="# of threads")
	parser.add_argument('--readNum', type=int, help="# of reads in fastq file")
	args = parser.parse_args()
	return args

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
	try:
		main()
	except SystemExit:
		raise
	except:
		print_exc()
		sys.exit(-1)
