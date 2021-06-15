#!/bin/python
# By Anoushka Joglekar 02.2021

"""  Takes in output of the 11-mer barcode-UMI matching algorithm
Output contains ordered list of 28mers with atleast 1 (or 2?) matching 11mers
Between ONT and 10X
Calculates pairwise edit distances on either strand per read
Outputs a file with lowest match edit distance on fwd and rev strand per read
"""

import os
import sys
import argparse
import pandas as pd
import Levenshtein as lv
from traceback import print_exc
import gzip

def perLine(line):
	""" Input can have 'NONE' on forward or reverse strand or both. Currently keeping
	cases where both have none while changing the distance to NA instead of zero.
	"""
	kmer_strings = [line[ix].split(',') for ix in [3,4,6,7]]
	read = line[0]
	gene = line[1]
	if 'NONE' not in kmer_strings[0]:
		lv_fwd = [lv.distance(kmer_strings[0][i],kmer_strings[1][i]) for i in range(len(kmer_strings[0]))]
		mF = min(lv_fwd)
		fwd_kmer = [kmer_strings[ix][lv_fwd.index(mF)] for ix in [0,1]]
	else:
		fwd_kmer = [kmer_strings[ix][0] for ix in [0,1]]
		mF = 'NA'
	if 'NONE' not in kmer_strings[2]:
		lv_rev = [lv.distance(kmer_strings[2][i],kmer_strings[3][i]) for i in range(len(kmer_strings[2]))]
		mR = min(lv_rev)
		rev_kmer = [kmer_strings[ix][lv_rev.index(mR)] for ix in [2,3]]
	else:
		rev_kmer = [kmer_strings[ix][0] for ix in [2,3]]
		mR = 'NA'
	fullOut = [read,gene,mF,fwd_kmer,mR,rev_kmer]
	""" fullOut has all Info on both, forward and reverse strands for further comparisons
	impInfo only has the minimum distance and k-mer. Cases where both have 'NONE' will be
	excluded later
	"""
	if mF != "NA" and mR != "NA":
		if mF < mR:
			out = [kmer_strings[ix][lv_fwd.index(min(lv_fwd))] for ix in [0,1]]
			d = lv.distance(out[0][0:16],out[1][0:16])
			impOut = [read, gene, mF, out[0], out[1], d]
		else:
			out = [kmer_strings[ix][lv_rev.index(min(lv_rev))] for ix in [2,3]]
			d = lv.distance(out[0][0:16],out[1][0:16])
			impOut = [read, gene, mR, out[0], out[1], d]
	elif mF == "NA" and mR == "NA":
		impOut = [read, gene, mR, 'NONE', 'NONE']
	elif mF == "NA" and mR != "NA":
		out = [kmer_strings[ix][lv_rev.index(mR)] for ix in [2,3]]
		d = lv.distance(out[0][0:16],out[1][0:16])
		impOut = [read, gene, mR, out[0], out[1], d]
	elif mF != "NA" and mR == "NA":
		out = [kmer_strings[ix][lv_fwd.index(mF)] for ix in [0,1]]
		d = lv.distance(out[0][0:16],out[1][0:16])
		impOut = [read, gene, mF, out[0], out[1], d]
	return([fullOut, impOut])


def main():
	""" Read in file and calculate pairwise edit dist.
	Kept chunking possible by having data be calculated in an
	indexed fashion. Seems unneccessary atm because v fast,
	i.e. 0.22 seconds per 10000 reads
	"""
	args = parse_args()
	readList = [a.decode('utf8').strip('\n').split('\t') for a in gzip.open(args.inFile,'rb')]
	""" can change range to chunks and add a diff function later
	"""
	multIter = [perLine(readList[ix]) for ix in range(len(readList))]
	data = [multIter[i][0] for i in range(len(multIter))]
	impData = [multIter[i][1] for i in range(len(multIter))]
	df = pd.DataFrame(data, columns = ['Read', 'Gene', 'minFwd', 'fwd_kmer', 'minRev', 'rev_kmer'])
	df.to_csv(args.outName+"_full.csv",sep="\t", index = False)
	os.system("cat %s | tr -d \"[]\' \" > tmp; mv tmp %s" % (args.outName+"_full.csv",args.outName+"_full.csv"))
	df_imp = pd.DataFrame(impData, columns = ['Read', 'Gene', 'minDist', '10x_kmer', 'ONT_kmer', 'bc_editDist'])
	df_imp = df_imp[df_imp.minDist!='NA']
	df_imp.to_csv(args.outName+"_min.csv",sep="\t", index = False)
	return()

def parse_args():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--inFile', type=str, help='Tab separated file from 11mer matching')
	parser.add_argument('--outName', type=str, help='output name')
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
