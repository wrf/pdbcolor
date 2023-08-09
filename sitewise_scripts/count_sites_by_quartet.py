#!/usr/bin/env python
#
# count_sites_by_quartet.py  created 2018-09-21

'''count_sites_by_quartet.py  last modified 2018-09-21
    count amino acid distributions around each branch of a phylogenetic quartet

count_sites_by_quartet.py -a matrix1.aln -l RAxML_perSiteLLs.matrix1.tab -q quartets.tab

    matrix can be in alternate formats (use -f), and gzipped

    generate the tabular sitewise results using:
sitewise_ll_to_columns.py RAxML_perSiteLLs.matrix1 > RAxML_perSiteLLs.matrix1.tab

    tabular likelihood output (-l) is assumed to look like:
site    T1      T2 ...
1       -12.345 -12.456 ...

    quartet table is tab delimited, for each species, indicated as A,B,C, or O
Human   A
Elephant    B
'''

import sys
import argparse
import time
import gzip
from collections import Counter,defaultdict
from Bio import AlignIO

def count_sites_for_quartets(fullalignment, alignformat, quartetspecies, toptreedict, lnldict, strongcutoff, skipgaps=True):
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

	# for each species, get index in the alignment
	quartetindex = defaultdict(list) # key is letter, value is list of indices
	for i, seqrec in enumerate(alignedseqs):
		species = str(seqrec.id)
		branch = quartetspecies[species]
		quartetindex[branch].append(i)
	if len(quartetindex) > 4:
		print >> sys.stderr, "WARNING: quartet has {} branches, should have only 4, check {}".format(len(quartetindex))

	strongcounter = 0
	print >> sys.stderr, "# determining AA for quartets", time.asctime()
	for i in range(al_length):
		aa_by_group = defaultdict( lambda: defaultdict(int) )
		for group,splist in quartetindex.iteritems():
			for sp in splist:
				aa_site = alignedseqs[sp,i]
				if skipgaps and ( aa_site=="-" or aa_site=="X" ):
					continue
				aa_by_group[group][aa_site] += 1
		dlnl = lnldict.get(i+1,0) # position are in alignment index, must offset by 1

		most_freq_O = sorted(aa_by_group["O"].items(), key=lambda x: x[1], reverse=True)
		most_freq_A = sorted(aa_by_group["A"].items(), key=lambda x: x[1], reverse=True)
		most_freq_B = sorted(aa_by_group["B"].items(), key=lambda x: x[1], reverse=True)
		most_freq_C = sorted(aa_by_group["C"].items(), key=lambda x: x[1], reverse=True)
		if most_freq_O and most_freq_A and most_freq_B and most_freq_C:

			if most_freq_O[0][0] == most_freq_A[0][0] and most_freq_B[0][0] == most_freq_C[0][0] and most_freq_O[0][0] != most_freq_C[0][0]:
				print >> sys.stdout, "{}\t{}\t{}\t{}\t{}\t{}\t{}".format( i+1, dlnl, toptreedict[i+1], aa_by_group["O"].items(), aa_by_group["A"].items(), aa_by_group["B"].items(), aa_by_group["C"].items() )
	#	if dlnl >= strongcutoff:
	#		strongcounter += 1
	#		print >> sys.stdout, "{}\t{}\t{}\t{}\t{}\t{}\t{}".format( i+1, dlnl, toptreedict[i+1], aa_by_group["O"].items(), aa_by_group["A"].items(), aa_by_group["B"].items(), aa_by_group["C"].items() )
		if i>=100000:
			break
	#	alignment_column = alignedseqs[:,i] # all letters per site
	#	numgaps = alignment_column.count("-")
	#	nogap_alignment_column = alignment_column.replace("-","").replace("X","") # excluding gaps
	#	num_nogap_taxa = len(nogap_alignment_column)
	#	aa_counter = Counter( nogap_alignment_column )
	#	mostcommonaa = aa_counter.most_common(1)[0][0]

	print >> sys.stderr, "# Counted {} strong sites".format( strongcounter )
	# no return

def read_tabular_3ln(lntabular):
	'''read tabular log-likelihood results and return a dict where keys are position and value is three-tree dlnL'''
	lnldict = {}
	toptreedict = {} # key is position, value is index of top tree 1 2 or 3
	linecounter = 0
	print >> sys.stderr, "# Reading log-likelihood by site from {}".format(lntabular), time.asctime()
	for line in open(lntabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			linecounter += 1
			if linecounter < 2:
				continue
			lsplits = line.split('\t')
			pos = int(lsplits[0])
			if lsplits[1]=="const": # encode constants as x
				lnldict[pos] = 0.0
				toptreedict[pos] = -1
			else: # for all other sites, do the calculation
				lnlfloats = map(float, lsplits[1:])
				toptree = lnlfloats.index(max(lnlfloats))
				rankedlnl = sorted(lnlfloats, reverse=True) # go from max likelihood to least
				dlnl = rankedlnl[0] - rankedlnl[1] # should always be greater than 0
				lnldict[pos] = dlnl
				toptreedict[pos] = toptree
	print >> sys.stderr, "# Found log-likelihood for {} sites".format( len(lnldict) ), time.asctime()
	return lnldict, toptreedict

def read_quartets(quartetfile):
	'''read quartet file and return a dict where key is species and value OABC group'''
	quartetdict = {}
	print >> sys.stderr, "# reading species for quartet from {}".format(quartetfile)
	for line in open(quartetfile,'r'):
		if line and line[0]!="#":
			species, branch = line.strip().split("\t")
			quartetdict[species] = branch
	return quartetdict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-l','--log-likelihood', help="tabular log-likelihood data file")
	parser.add_argument('-q','--quartet', help="tabular file of species to quartet branch")
	parser.add_argument('-t','--strong', type=float, default=0.5, help="cutoff threshold for strong sites in RAxML or phylobayes [0.5]")
	args = parser.parse_args(argv)

	species_to_quartet = read_quartets(args.quartet)
	lnldict, toptrees = read_tabular_3ln(args.log_likelihood)
	count_sites_for_quartets(args.alignment, args.format, species_to_quartet, toptrees, lnldict, args.strong)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
