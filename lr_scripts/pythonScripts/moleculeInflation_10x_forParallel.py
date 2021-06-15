#! /bin/python3

# By Anoushka Joglekar 01.25.2021
""" For ONT data we want to collapse UMIs and avoid
situations where an error converts one UMI to another
and inflates molecule counts """

import sys
import argparse
import pandas as pd
import Levenshtein as lv
import contextlib
import joblib
from tqdm import tqdm
from joblib import Parallel, delayed
from traceback import print_exc
import time
from datetime import timedelta
import re

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
	"""Context manager to patch joblib to report into tqdm progress bar given as argument"""
	class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
		def __init__(self, *args, **kwargs):
			super().__init__(*args, **kwargs)

		def __call__(self, *args, **kwargs):
			tqdm_object.update(n=self.batch_size)
			return super().__call__(*args, **kwargs)
	old_batch_callback = joblib.parallel.BatchCompletionCallBack
	joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
	try:
		yield tqdm_object
	finally:
		joblib.parallel.BatchCompletionCallBack = old_batch_callback
		tqdm_object.close()


def processPerBarcode(bc, gene, umiDict_ont, umiDict_tx):
	""" For each barcode of a gene in the ONT data,
	get corresponding UMIs from both platforms.
	If there is no barcode in PB, we skip, else we
	find the minimum Levenshtein distance
	"""
	ont = umiDict_ont['UMI'][bc].split(',')
	try:
		tx = umiDict_tx['UMI'][bc].split(',')
		l = []
		for umi in ont:
			if umi in tx:
				l.append([umi,0,umi])
			else:
				ld = [lv.distance(umi,u) for u in tx]
				m = min(ld)
				l.append([umi,m,tx[ld.index(m)]])
	except KeyError:
		l = [[umi,'NA',umi] for umi in ont]
	l2 = []
	for umi in ont:
		if len(ont) > 1:
#			ont_prime = ont.copy()
			ont_prime = ont
			ont_prime.remove(umi)
			ld = [lv.distance(umi,u) for u in ont_prime if umi != u]
			m = min(ld)
			l2.append([umi,m,ont_prime[ld.index(m)]])
		else:
			l2.append([umi,'NA',umi])
	a = pd.DataFrame(l, columns = ['umi','dist','alt_umi'])
	b = pd.DataFrame(l2, columns = ['umi','self_dist','self_umi'])
	bc_df = pd.merge(a,b)
	bc_df['Gene'] = gene
	bc_df['Barcode'] = bc
	return(bc_df)

def subsetDF_toDict(df, gene):
	""" Subset a dataframe by gene and convert to
	dictionary by using the barcodes as keys. The
	dicts are sorted in a descending order of count, but
	that's functionality for later """
	pG = df[df['Gene']==gene].\
		sort_values(by = ['count'], ascending = False).\
		reset_index(drop = True)
	umiDict = pG.drop('Gene', axis = 1).\
		groupby('Barcode').\
		agg(lambda x: ','.join(x)).to_dict()
	return(umiDict)

def getClosestUMI(gene, ont_df, txdf):
	""" Iterate the barcode processing over
	all barcodes of a specific gene """
	umiDict_tx = subsetDF_toDict(txdf,gene)
	umiDict_ont = subsetDF_toDict(ont_df,gene)
	barcodeList = list(umiDict_ont['UMI'].keys())
	bcDF_list = []
	for barcode in barcodeList:
		bcDF_list.append(processPerBarcode(barcode, gene, umiDict_ont, umiDict_tx))
	geneDF = pd.concat(bcDF_list).reset_index(drop = True)
	return(geneDF)

def processAI(inFile,cols,ontFile):
	""" Read in AllInfo files (could be complete or
	incomplete reads) and convert to a DF b.c. that's
	easier to work with """
	if ontFile == "T":
		readList = [[re.split('\t|\.',a.strip('\n'))[i] for i in cols] for a in open(inFile,'r')]
	else:
		readList = [[re.split('\t',a.strip('\n'))[i] for i in cols] for a in open(inFile,'r')]
	df = pd.DataFrame(readList, columns = ['Gene','Barcode','UMI'])
	df['count'] = df.groupby(['Gene','Barcode']).transform('count')
	cond_df = df.groupby(['Gene','Barcode','UMI']).agg(['count']).reset_index(drop=False)
	cond_df.columns = ['Gene','Barcode','UMI','count']
	return(cond_df)

def main():
	""" Read in input, iterate over each gene in a parallel
	fashion
	"""
	args = parse_args()
	t0 = time.time()
	txdf = processAI(args.tx,[1,3,4],'F')
	ont_df = processAI(args.ont,[1,4,5],'T')
	t1 = time.time()
	print("Processed input in (hh:mm:ss.ms) ",str(timedelta(seconds=t1-t0)))
	geneList = list(set([gene for gene in ont_df['Gene']]))
	g_df = []
	if "_00" in args.ont:
		with tqdm_joblib(tqdm(desc="# genes done", total=len(geneList))) as progress_bar:
			nested = Parallel(n_jobs=args.threads,prefer="threads")(delayed(getClosestUMI)(gene, ont_df, txdf) for gene in geneList)
	else:
			nested = Parallel(n_jobs=args.threads,prefer="threads")(delayed(getClosestUMI)(gene, ont_df, txdf) for gene in geneList)
	fullDF = pd.concat(nested)
	fullDF.to_csv(args.outName+".csv",sep="\t", index = False)
	return()

def parse_args():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--ont', type=str, help='AllInfo file from ONT data')
	parser.add_argument('--outName', type=str, help='output name')
	parser.add_argument('--tx', type=str, help='10x whitelist')
	parser.add_argument('--threads', type=int, help='# of threads to use')
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


