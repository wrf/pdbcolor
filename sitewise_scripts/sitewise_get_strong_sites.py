#!/usr/bin/env python
#
# sitewise_get_strong_sites.py  created 2018-02-01

'''sitewise_get_strong_sites.py  last modified 2018-03-15
    subset a supermatrix using only strong-sites identified by site-wise RAxML

sitewise_get_strong_sites.py -a matrix1.aln -l RAxML_perSiteLLs.matrix1.tab -o strong_alignment.aln

    matrix can be in alternate formats (use -f), and gzipped

    for phylobayes sitelogl values, it is advisable to set -m to 1.0

    generate the tabular sitewise results using:
sitewise_ll_to_columns.py RAxML_perSiteLLs.matrix1 > RAxML_perSiteLLs.matrix1.tab
'''

import sys
import argparse
import time
import gzip
from Bio import AlignIO
from numpy import median

def make_strong_alignment(fullalignment, alignformat, valsbysite, mindlnl, makeweak):
	'''read large alignment, and return a new alignment of only strong sites'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	numtaxa = len(alignedseqs)
	allength = alignedseqs.get_alignment_length()

	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( numtaxa, allength )
	strongcounter = 0
	nullcounter = 0
	newalign = alignedseqs[:,0:0] # start with blank alignment

	for i in range(allength): # sites begin at 0 while lnl begins at 1
		if i>0 and i%10000==0:
			sys.stderr.write(".")
		if makeweak:
			if valsbysite[i+1] >= mindlnl:
				strongcounter += 1
				continue
			nullcounter += 1
			newalign += alignedseqs[:,i:i+1]
		else:
			if valsbysite[i+1] >= mindlnl:
				strongcounter += 1
				newalign += alignedseqs[:,i:i+1]
	if nullcounter:
		print >> sys.stderr, "# New alignment contains {} weak sites, skipped {} strongs sites".format( nullcounter, strongcounter )
	elif strongcounter:
		print >> sys.stderr, "# New alignment contains {} strong sites".format( strongcounter )
	return newalign

def read_tabular_ln(lntabular, treelist):
	'''read tabular log-likelihood results and return a dict where keys are position and value is abs dlnL'''
	lndict = {}
	linecounter = 0
	print >> sys.stderr, "# Using columns {} and {} for topologies".format(*treelist)
	print >> sys.stderr, "# Reading log-likelihood by site from {}".format(lntabular), time.asctime()
	for line in open(lntabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			linecounter += 1
			if linecounter < 2:
				continue
			lsplits = line.split('\t')
			pos = int(lsplits[0]) # sites begin at 1
			lnlfloats = [float(n) for n in lsplits[1:]]
			dlnl = max(lnlfloats) - median(lnlfloats) # should always be greater than 0
			lndict[pos] = dlnl
	print >> sys.stderr, "# Found log-likelihood for {} sites".format( len(lndict) ), time.asctime()
	return lndict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-l','--log-likelihood', help="tabular log-likelihood data file from RAxML")
	parser.add_argument('-m','--dlnl-minimum', default=0.5, type=float, help="minimum difference in lnL [0.5]")
	parser.add_argument('-o','--output', help="name of output file", required=True)
	parser.add_argument('-t','--trees', default="1,2", help="two columns of trees 1 and 2, as comma-separated ints [1,2]")
	parser.add_argument('-w','--weak', action="store_true", help="get weak sites, meaning exclude strong sites")
	args = parser.parse_args(argv)

	treelist = [int(i) for i in args.trees.split(",")]
	valsbysite = read_tabular_ln(args.log_likelihood, treelist)
	strongalignment = make_strong_alignment(args.alignment, args.format, valsbysite, args.dlnl_minimum, args.weak)
	AlignIO.write(strongalignment, args.output, args.format)
	print >> sys.stderr, "# New alignment written to {}".format(args.output), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
