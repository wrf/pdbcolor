#!/usr/bin/env python
#
# pdb_log_likelihood.py v1 2018-02-01

'''pdb_log_likelihood.py  last modified 2018-02-13

pdb_log_likelihood.py -a PRP4B_HUMAN.aln -p 4ian.pdb -s PRP4B_HUMAN > 4ian_w_lnl.pdb

    for PDB files that contain multiple proteins (with fewer than 100k atoms)
    multiple alignments and seq IDs can be given by listing after -a and -s
    the order of the files must match the sequence IDs

pdb_log_likelihood.py -p 2o8b.pdb -a 10543-11140-MSH2_HUMAN.aln 51199-52049-MSH6_HUMAN.aln -s MSH2_HUMAN MSH6_HUMAN > 2o8b_w_both_lnl.pdb

within a PDB file, fields for atoms are:
Record name      Residue         Position as X Y Z
          Atom serial number                            Occupancy
             Atom    Chain                                   temperatureFactor
                        Residue number
such as:
ATOM      1  N   PRO A  46       0.739  40.031  44.896  1.00 0.842           N

  here, the temperatureFactor will be replaced with the delta-log-likelihood
  as calculated by:

blast_to_align_pairs.py -b simion_taxa/hsapiens_vs_uniprot_blastp.tab -q simion_taxa/Homo_sapiens.fasta.nogaps -s ~/db/human_uniprot.fasta -r simion_taxa/Homo_sapiens.fasta -l RAxML_perSiteLLs.simion2017.tab

residue number must match the alignment, not the position in the model
meaning even if the PDB file starts with residue 20, if the first 19 were
disordered or cleaved, the sequence still must start with residue 1
'''

import sys
import time
import argparse
from collections import Counter
from Bio import AlignIO

def get_alignment_values(alignmentlist, alignformat, targetidlist):
	scoreindex_dict = {} # dict of dicts, where key is seqID, value is dict where key is position
	for alignment,target_seqid in zip(alignmentlist,targetidlist):
		print >> sys.stderr, "# Reading alignment from {}".format( alignment )
		alignment = AlignIO.read( alignment, alignformat )

		al_length = alignment.get_alignment_length()
		num_taxa = len(alignment)

		print >> sys.stderr, "# Alignment contains {} sequences for {} sites, including gaps".format( num_taxa, al_length )

		targetseq = None
		# scores can be likelihood or heteropecilly
		sitescores = None
		scoretype = None
		for seqrec in alignment:
			try: # this is for Uniprot IDs
				shortid = seqrec.id.split("|")[2]
			except IndexError:
				shortid = None
			if seqrec.id==target_seqid or shortid==target_seqid:
				targetseq = seqrec.seq
			elif seqrec.id=="Likelihood_score" or seqrec.id=="Heteropecilly_score":
				sitescores = seqrec.seq
				scoretype = seqrec.id
		if targetseq is None:
			print >> sys.stderr, "# ERROR: CANNOT FIND SEQUENCE {}, CHECK OPTION -s OR ALIGNMENT".format( target_seqid )
			return None
		if sitescores is None:
			print >> sys.stderr, "# ERROR: CANNOT FIND SITEWISE SCORES {}, CHECK ALIGNMENT".format( target_seqid )
			return None

		index_to_score = {}
		nongapcount = 0
		targetcount = 0 # keep track of position in target sequence for all non-gap letters
		if scoretype=="Likelihood_score":
			for i in range(al_length):
				targetletter = targetseq[i]
				lnlvalue = sitescores[i]
				if targetletter != "-": # meaning anything except gaps in target, though there should not be any as it is the reference sequence
					targetcount += 1
					if lnlvalue == "-":
						rankedscore = -1
					elif lnlvalue == "x" or lnlvalue == "X":
						nongapcount += 1
						rankedscore = 16
					else:
						nongapcount += 1
						rankedscore = int(lnlvalue,16)
					index_to_score[targetcount] = rankedscore
			print >> sys.stderr, "# Found likelihood for {} sites for {}".format( nongapcount, target_seqid )
			scoreindex_dict[target_seqid] = index_to_score
		elif scoretype=="Heteropecilly_score":
			for i in range(al_length):
				targetletter = targetseq[i]
				hpvalue = sitescores[i]
				if targetletter != "-":
					targetcount += 1
					if hpvalue == "-":
						rankedscore = -1
					elif hpvalue == "c":
						nongapcount += 1
						rankedscore = 10
					elif hpvalue == "C":
						nongapcount += 1
						rankedscore = 11
					else:
						nongapcount += 1
						rankedscore = int(hpvalue)
					index_to_score[targetcount] = rankedscore
			print >> sys.stderr, "# Found heteropecilly for {} sites for {}".format( nongapcount, target_seqid )
			scoreindex_dict[target_seqid] = index_to_score
	return scoreindex_dict

def rewrite_pdb(pdbfile, seqidlist, scoredict, wayout, forcerecode, colorgaps, heterocolors):
	print >> sys.stderr, "# Reading PDB from {}".format(pdbfile)
	atomcounter = 0
	residuecounter = {}
	keepchains = {} # dict where key is chain and value is seqid
	defaultchain = True # flag for whether DBREF occurs at all
	for line in open(pdbfile,'r'):
		# records include:
		# HEADER TITLE COMPND SOURCE AUTHOR REVDAT JRNL REMARK
		# DBREF SEQRES HET HETNAM FORMUL HELIX SHEET SSBOND LINK CISPEP SITE ATOM CONECT

		#COLUMNS        DATA  TYPE    FIELD        DEFINITION
		#-------------------------------------------------------------------------------------
		# 1 -  6        Record name   "ATOM  "
		# 7 - 11        Integer       serial       Atom  serial number.
		#13 - 16        Atom          name         Atom name.
		#17             Character     altLoc       Alternate location indicator.
		#18 - 20        Residue name  resName      Residue name.
		#22             Character     chainID      Chain identifier.
		#23 - 26        Integer       resSeq       Residue sequence number.
		#27             AChar         iCode        Code for insertion of residues.
		#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		#55 - 60        Real(6.2)     occupancy    Occupancy.
		#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		#77 - 78        LString(2)    element      Element symbol, right-justified.
		#79 - 80        LString(2)    charge       Charge  on the atom.

		record = line[0:6].strip()
		if record=="DBREF":
			defaultchain = False
			proteinid = line[42:56].strip()
			for seqid in seqidlist:
				if seqid.find(proteinid)>-1:
					chaintarget = line[12]
					print >> sys.stderr, "### keeping chain {} for sequence {}".format( chaintarget, proteinid )
					keepchains[chaintarget] = proteinid
		# for all other lines, check for ATOM or not
		if record=="ATOM": # skip all other records
			chain = line[21]
			residue = int( line[22:26] )
			if defaultchain or forcerecode or chain in keepchains:
				atomcounter += 1
				if defaultchain or forcerecode: # assume only one seqid
					score = scoredict[seqidlist[0]].get(residue,0.00)
				else:
					score = scoredict[keepchains[chain]].get(residue,0.00)
				if score:
					residuecounter[residue] = True
			else: # meaning in another chain, so color as insufficient
				score = 99
			newline = "{}{:6.2f}{}".format( line[:60], score, line[66:].rstrip() )
			print >> wayout, newline
		else:
			print >> wayout, line.strip()
	print >> sys.stderr, "# Recoded values for {} atoms in {} residues".format(atomcounter, len(residuecounter) )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", nargs="*", help="multiple sequence alignment", required=True)
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	parser.add_argument("-s","--sequence", nargs="*", help="sequence ID for PDB", required=True)
	parser.add_argument("--color-gaps", action="store_true", help="use dark colors to uniquely color gaps of up to 6 chains")
	parser.add_argument("--force-recode", action="store_true", help="force recoding regardless of chain")
	parser.add_argument("--heteroatoms", action="store_true", help="color heteroatoms")
	args = parser.parse_args(argv)

	if len(args.alignment) != len(args.sequence):
		print >> sys.stderr, "ERROR: {} ALIGNMENTS FOR {} SEQUENCES, MUST BE EQUAL, CHECK -a AND -s".format(len(args.alignment), len(args.sequence)), time.asctime()

	scoredict = get_alignment_values( args.alignment, args.format, args.sequence)
	if scoredict: # indicating that the sequence was found and something was calculated
		rewrite_pdb(args.pdb, args.sequence, scoredict, wayout, args.force_recode, args.color_gaps, args.heteroatoms)
	else:
		sys.exit("# CANNOT CALCULATE LIKELIHOODS, EXITING")

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
