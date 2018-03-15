#!/usr/bin/env python
#
# sitewise_to_genewise_lnl.py  created 2018-03-01

'''sitewise_to_genewise_lnl.py  last modified 2018-03-01
    recode constant sites as null from site-wise RAxML tabular output

sitewise_to_genewise_lnl.py -l RAxML_perSiteLLs.matrix1.tab -p partitions.txt > matrix1_per_gene_lnl.tab

    generate the tabular sitewise results using:
sitewise_to_genewise_lnl.py RAxML_perSiteLLs.matrix1 > RAxML_perSiteLLs.matrix1.tab
'''

import sys
import argparse
import time
import gzip
from collections import Counter,defaultdict
from Bio import AlignIO

def get_partitions(partitionfile):
	'''read comma-delimited partition information and return a list of tuples'''
	partitions = [] # list of tuples of intervals
	for line in open(partitionfile,'r'):
		line = line.strip()
		if line:
			blocks = line.split(",") # split "1:136,137:301,..." into ['1:136', '137:301',...]
			for block in blocks:
				alignindex = tuple( int(i) for i in block.split(":") ) # split '1:136' into ( 1,136 )
				partitions.append(alignindex)
	print >> sys.stderr, "# read {} partitions from {}".format(len(partitions), partitionfile), time.asctime()
	return partitions

def sum_sites_by_gene(lntabular, partitions, wayout):
	'''read tabular log-likelihood results and sum by gene partitions'''

	genesums = {} # keys are partition tuples, values are list of lnL sums

	# cycle through the parts, create dict for each site to its partition
	sitetopart = {} # key is site, value is partition
	for part in partitions:
		for i in range(part[0], part[1]+1):
			sitetopart[i] = part

	linecounter = 0
	print >> sys.stderr, "# Reading log-likelihood by site from {}".format(lntabular), time.asctime()
	for line in open(lntabular,'r'):
		line = line.strip()
		if line: # ignore blank lines
			linecounter += 1
			lsplits = line.split('\t')
			if linecounter < 2:
				numtrees = len(lsplits)-1
				for part in partitions:
					genesums[part] = [0.0] * numtrees
				print >> wayout, line
				continue
			pos = int(lsplits[0]) # sites begin at 1
			targetpart = sitetopart[pos]
			for i,swl in enumerate(lsplits[1:]):
				if swl=="const":
					genesums[targetpart][i] += 0.0
				else:
					genesums[targetpart][i] += float(swl)
	print >> sys.stderr, "# Read {} lines".format( linecounter ), time.asctime()

	for part in partitions:
		print >> wayout, "{}\t{}".format("{}-{}".format(*part), "\t".join([str(x) for x in genesums[part]]))

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-l','--log-likelihood', help="tabular log-likelihood data file from RAxML")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	args = parser.parse_args(argv)

	partitions = get_partitions(args.partition)
	sum_sites_by_gene(args.log_likelihood, partitions, wayout)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
