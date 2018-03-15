#!/usr/bin/env python
#
# sitewise_recode_constant_sites.py  created 2018-02-02

'''sitewise_recode_constant_sites.py  last modified 2018-03-12
    recode constant sites as null from site-wise RAxML or phylobayes tabular output

sitewise_recode_constant_sites.py -a matrix1.aln -l RAxML_perSiteLLs.matrix1.tab > RAxML_perSiteLLs.matrix1_w_const.tab

    matrix can be in alternate formats (use -f), and gzipped

    generate the tabular sitewise results using:
sitewise_ll_to_columns.py RAxML_perSiteLLs.matrix1 > RAxML_perSiteLLs.matrix1.tab

    tabular likelihood output (-l) is assumed to look like:
site    T1      T2 ...
1       -12.345 -12.456 ...
'''

import sys
import argparse
import time
import gzip
from collections import Counter
from Bio import AlignIO

def determine_constant_sites(fullalignment, alignformat):
	'''read large alignment, and return a dict where keys are constant positions'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	numtaxa = len(alignedseqs)
	al_length = alignedseqs.get_alignment_length()
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( numtaxa, al_length )

	constsites = {}
	print >> sys.stderr, "# determining constant sites", time.asctime()
	for i in range(al_length):
		alignment_column = alignedseqs[:,i] # all letters per site
		numgaps = alignment_column.count("-")
		nogap_alignment_column = alignment_column.replace("-","").replace("X","") # excluding gaps
		num_nogap_taxa = len(nogap_alignment_column)
		aa_counter = Counter( nogap_alignment_column )
		mostcommonaa = aa_counter.most_common(1)[0][0]

		if len(aa_counter)==1: # meaning site has more than 1 possible AA, so use lnL
			constsites[i+1] = True # use index plus 1 to match positions
			mostcommoncount = aa_counter.most_common(1)[0][1]
	print >> sys.stderr, "# Counted {} constant sites".format( len(constsites) )
	return constsites

def read_tabular_ln(lntabular, constsites, wayout):
	'''read tabular log-likelihood results and print a modified table where constants sites are recoded'''
	linecounter = 0
	recodecount = 0
	print >> sys.stderr, "# Reading log-likelihood by site from {}".format(lntabular), time.asctime()
	for line in open(lntabular,'r'):
		line = line.strip()
		if line: # ignore blank lines
			linecounter += 1
			if linecounter < 2:
				print >> wayout, line
				continue
			lsplits = line.split('\t')
			pos = int(lsplits[0]) # sites begin at 1
			if constsites.get(pos, False): # meaning site is constant, so recode
				recodecount += 1
				recodedsites = ["const"] * len(lsplits[1:])
				print >> wayout, "{}\t{}".format( pos, "\t".join(recodedsites) )
			else:
				print >> wayout, line
	print >> sys.stderr, "# Read {} lines, recoded {} sites".format( linecounter,recodecount ), time.asctime()

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-l','--log-likelihood', help="tabular log-likelihood data file from RAxML")
	args = parser.parse_args(argv)

	constsitedict = determine_constant_sites(args.alignment, args.format)
	read_tabular_ln(args.log_likelihood, constsitedict, wayout)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
