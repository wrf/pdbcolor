#!/usr/bin/env python
#
# blast_to_align_pairs.py
# pdbcolor version v2018-09-21
# v1.1 2023-02-01 python3 update
# v1.2 2025-12-27 indicate likelihood3 for 3-tree calculations

'''blast_to_align_pairs.py  last modified 2025-12-27
    this is the script that generates the many alignments of each protein
    where one line is the alignment from the supermatrix,
    one line is the original protein (untrimmed),
    and one line is the log-likelihood or heteropecilly matched with the aligned portion

>Homo_sapiens_14885-15764
-GFLKLIEIENFKSYKGRQIIGPFQ-FTAIIGPNGSGKSNLMDAISFVLGEKTSNLRVK-
>sp|Q14683|SMC1A_HUMAN Structural maintenance of chromosomes protein 1A OS=Homo sapiens GN=SMC1A PE=1 SV=2
MGFLKLIEIENFKSYKGRQIIGPFQRFTAIIGPNGSGKSNLMDAISFVLGEKTSNLRVKT
>Likelihood3_score
-66x00000000x0x006060xx06-x3030x0xx0xxxx00x000x000636660x66-

blast_to_align_pairs.py -q query_prots.fasta -s prot_db.fasta -b blastp.tab -d pair_dir -r Homo_sapiens.fasta

    make tabular blast output from:

blastp -query query_prots.fasta -db prot_db.fasta -outfmt 6 -max_target_seqs 1 > blastp.tab

    alignment reference proteins (still aligned) can be extracted with:

split_supermatrix_to_taxa.py -a matrix.phy -p partitions.txt -d taxa -f phylip-relaxed
'''

import sys
import os
import argparse
import time
import subprocess
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq

def read_tabular_hp(hptabular):
	'''read tabular heteropecilly results and return a dict where keys are position and value is HP'''
	hpdict = {}
	print( "# Reading heteropecilly by site from {}  {}".format(hptabular, time.asctime() ), file=sys.stderr )
	for line in open(hptabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			pos, hp = line.split('\t')[0:2] # only take first two columns
			pos = int(pos)
			if hp == "10":
				hp = "c" # redefine constant sites as C instead of 10 for string
			elif hp == "11":
				hp = "C"
			hpdict[pos] = hp
	print( "# Found heteropecilly for {} sites  {}".format( len(hpdict) , time.asctime() ), file=sys.stderr )
	return hpdict

def read_tabular_3ln(lntabular, strongsitecutoff):
	'''read tabular log-likelihood results and return a dict where keys are position and value is three-tree dlnL'''
	lnldict = {}
	linecounter = 0
	print( "# Reading log-likelihood by site from {}  {}".format(lntabular, time.asctime() ), file=sys.stderr )
	for line in open(lntabular,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore blank and comment lines
			linecounter += 1
			lsplits = line.split('\t')
			if linecounter < 2:
				if len(lsplits) < 4:
					print( "# WARNING: file has {} columns, not 4".format( len(lsplits) ), file=sys.stderr )
				continue
			pos = int(lsplits[0])
			if lsplits[1]=="const": # encode constants as x
				lnldict[pos] = "x"
			else: # for all other sites, do the calculation
				lnlfloats = list(map(float, lsplits[1:]))
				# top tree for each site before values are sorted
				toptree = lnlfloats.index(max(lnlfloats))
				# regardless of how many trees, sort them from max likelihood to least
				rankedlnl = sorted(lnlfloats, reverse=True)
				dlnl = rankedlnl[0] - rankedlnl[1] # should always be greater than 0
				# was calculated before as mean of T1-T2, T2-T3, T1-T3
				# but should instead be max likelihood - second best alternative
				# to avoid cases where min is much worse than max and med
				if dlnl < strongsitecutoff: # should be 0.5 by default
					adjdlnl = 0
				else: # meaning above noise threshold
					if dlnl < 2.0: # strong site, not very strong
						adjdlnl = 1
					else: # very strong sites
						adjdlnl = 2
				adjdlnl = adjdlnl + (3 * toptree) # assumes 3 trees, though would work with 2
				if adjdlnl >= 9:
					adjdlnl = 9
				elif adjdlnl <= 0:
					adjdlnl = 0
				lnldict[pos] = str(adjdlnl)

	print( "# Found log-likelihood for {} sites".format( len(lnldict) , time.asctime() ), file=sys.stderr )
	return lnldict

def run_mafft(MAFFT, rawseqsfile):
	'''generate multiple sequence alignment from fasta and return MSA filename'''
	aln_output = "{}.aln".format(os.path.splitext(rawseqsfile)[0] )
	aligner_args = [MAFFT, "--auto", "--quiet", rawseqsfile]
	print( "{}\n{}".format(time.asctime(), " ".join(aligner_args) ), file=sys.stderr )
	with open(aln_output, 'w') as msa:
		subprocess.call(aligner_args, stdout=msa)
	print( "# alignment of {} completed  {}".format(aln_output, time.asctime() ), file=sys.stderr )
	if os.path.isfile(aln_output):
		return aln_output
	else:
		raise OSError("Cannot find expected output file {}".format(aln_output) )

def make_aligned_score_string(aln_file, scoredict, refprotdict):
	'''from partial alignment to full protein, return string of sitewise-lnl or heteropecilly for partial alignment'''
	print( "# Reading alignment from {}".format( aln_file ), file=sys.stderr )
	alignment = AlignIO.read( aln_file, "fasta" )
	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)
	print( "# Alignment contains {} taxa for {} sites, including gaps".format( num_taxa, al_length ), file=sys.stderr )
	rangeindex = [ int(x) for x in os.path.basename(aln_file).split("-")[0:2] ]

	# remove gapped sites from heteropecilly based on aligned and trimmed protein
	score_no_gaps = ""
	hpoffset = 0
	refid = alignment[0].id # should always be first seq in alignment
	for i,site in enumerate(str(refprotdict[refid].seq)):
		if site!="-": # if site is not gap, keep heteropecilly score
			adjindex = rangeindex[0] + i # + hpoffset
			score_no_gaps += scoredict[adjindex]
	#	else: # increase offset for each gap
		#	hpoffset += 1

	# assign heteropecilly score to columns without gaps in alignment
	aligned_scorestring = ""
	refoffset = 0
	for i in range(al_length):
		alignment_column = alignment[:,i] # all letters per site
		if alignment_column[0]=="-": # if first character is a gap, assign gap for heteropecilly
			aligned_scorestring += "-"
		else: # otherwise use value from score_no_gaps
			aligned_scorestring += score_no_gaps[refoffset]
			refoffset += 1
	return aligned_scorestring

def make_pairs_from_blast(blastfile, querydict, subjectdict, new_aln_dir, mafftbin, scoredict, refdict, fastascorestring ):
	'''iterate through blast hits, generate fasta files of each query subject pair and align them'''
	qspairs = {} # key is query, value is subject, assuming only one hit but multiple HSPs
	linecounter = 0
	print( "# parsing blast hits from {}  {}".format(blastfile, time.asctime() ), file=sys.stderr )
	for line in open(blastfile,'r'):
		line = line.strip()
		if line and line[0]!="#": # skip empty and comment lines
			linecounter += 1
			lsplits = line.split("\t")
			query = lsplits[0]
			subject = lsplits[1]
			qspairs[query] = subject
	print( "# counted {} blast hits for {} query proteins  {}".format( linecounter, len(qspairs) , time.asctime() ), file=sys.stderr )
	filecounter = 0
	for q,s in qspairs.items():
		fasta_filename = "{}-{}.fasta".format(q.split("_")[-1],s.split("|")[-1])
		fasta_output = os.path.join(new_aln_dir,fasta_filename)
		with open(fasta_output,'w') as fo:
			fo.write( querydict[q].format("fasta") )
			if str(subjectdict[s].seq).count("U"): # check for disallowed selenocysteine
				subjectdict[s].seq = Seq(str(subjectdict[s].seq).replace("U","C"))
			fo.write( subjectdict[s].format("fasta") )
		aln_file = run_mafft(mafftbin, fasta_output)
		if aln_file: # check if file was made
			filecounter += 1
			if refdict and scoredict: # if alignment proteins and heteropecilly are given
				hpstring = make_aligned_score_string(aln_file, scoredict, refdict)
				with open(aln_file, 'a') as af:
					print( ">{}\n{}".format( fastascorestring, hpstring ), file=af )
	print( "# generated {} alignments  {}".format(filecounter, time.asctime() ), file=sys.stderr )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b','--blast', help="tabular blastp results")
	parser.add_argument('-d','--directory', default="./blast_alignments", help="directory for new alignments [autonamed]")
	parser.add_argument('--mafft', default="mafft", help="path to mafft binary [default is in PATH]")
	parser.add_argument('-q','--query', help="fasta file of blast query proteins from supermatrix")
	parser.add_argument('-s','--subject', help="fasta file of blast database proteins; this was the -db file for blastp")
	parser.add_argument('-r','--reference', help="fasta file of proteins from the alignment with gaps, only needed if log-likelihood or heteropecilly is computed")
	parser.add_argument('-t','--strong', type=float, default=0.5, help="cutoff threshold for strong sites in RAxML or phylobayes [0.5]")
	parser.add_argument('-l','--log-likelihood', help="tabular three-tree log-likelihood data file from RAxML")
	parser.add_argument('-p','--heteropecilly', help="tabular heteropecilly data file")
	args = parser.parse_args(argv)

	### DIRECTORY FOR NEW OUTPUT
	new_aln_dir = os.path.abspath("{}_{}".format(args.directory, time.strftime("%Y%m%d-%H%M%S") ) )
	if not os.path.exists(new_aln_dir):
		os.mkdir(new_aln_dir)
		print( "# Creating directory {}  {}".format(new_aln_dir, time.asctime() ), file=sys.stderr )
	elif os.path.isdir(new_aln_dir):
		print( "# Using directory {}  {}".format(new_aln_dir, time.asctime() ), file=sys.stderr )

	print( "# Reading sequences from {}  {}".format(args.query, time.asctime() ), file=sys.stderr )
	querydict = SeqIO.to_dict(SeqIO.parse(args.query,"fasta"))
	print( "# Counted {} sequences  {}".format( len(querydict) , time.asctime() ), file=sys.stderr )
	print( "# Reading sequences from {}  {}".format(args.subject, time.asctime() ), file=sys.stderr )
	subjectdict = SeqIO.to_dict(SeqIO.parse(args.subject,"fasta"))
	print( "# Counted {} sequences  {}".format( len(subjectdict) , time.asctime() ), file=sys.stderr )
	if args.reference:
		print( "# Reading sequences from {}  {}".format(args.reference, time.asctime() ), file=sys.stderr )
		refdict = SeqIO.to_dict(SeqIO.parse(args.reference,"fasta"))
		print( "# Counted {} sequences  {}".format( len(refdict) , time.asctime() ), file=sys.stderr )
	else:
		refdict = None

	valsbysite = None # by default is None
	fastascorestring = None
	if args.heteropecilly:
		valsbysite = read_tabular_hp(args.heteropecilly, args.strong)
		fastascorestring = "Heteropecilly_score"
	elif args.log_likelihood:
		valsbysite = read_tabular_3ln(args.log_likelihood, args.strong)
		fastascorestring = "Likelihood3_score"

	make_pairs_from_blast(args.blast, querydict, subjectdict, new_aln_dir, args.mafft, valsbysite, refdict, fastascorestring)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
