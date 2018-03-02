#!/usr/bin/env python
#
# pdb_site_identity.py v1 2017-07-25

'''pdb_site_identity.py  last modified 2018-01-08

pdb_site_identity.py -a mox_all.aln -s DOPO_HUMAN -p 4zel.pdb > 4zel_w_scores.pdb

within a PDB file, fields for atoms are:
Record name      Residue         Position as X Y Z
          Atom serial number                            Occupancy
             Atom    Chain                                   temperatureFactor
                        Residue number
such as:
ATOM      1  N   PRO A  46       0.739  40.031  44.896  1.00 0.842           N

here, the temperatureFactor will be replaced with either the identity %
or an integer for the ranked identity score, where bins are:
0, 50, 60, 70, 80, 90, 95, 98, 100%
corresponding to scores of 1 to 9

residue number must match the alignment, not the position in the model
meaning even if the PDB file starts with residue 20, if the first 19 were
disordered or cleaved, the sequence still must start with residue 1
'''

import sys
import time
import argparse
from collections import Counter
from Bio import AlignIO

def get_alignment_identity(alignmentlist, alignformat, targetidlist, gapcutoff, store_identity=False):
	identindex_dict = {} # key is target seqid, value is dict of identities
	for alignment, target_seqid in zip(alignmentlist,targetidlist):
		print >> sys.stderr, "# Reading alignment from {}".format( alignment )
		alignment = AlignIO.read( alignment, alignformat )

		al_length = alignment.get_alignment_length()
		num_taxa = len(alignment)

		print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( num_taxa, al_length )

		targetseq = None
		for seqrec in alignment:
			if seqrec.id==target_seqid:
				targetseq = seqrec.seq
		if targetseq is None:
			print >> sys.stderr, "# ERROR: CANNOT FIND SEQUENCE {}, CHECK OPTION -s OR ALIGNMENT".format( target_seqid )
			return None

		index_to_identity = {}
		targetcount = 0 # keep track of position in target sequence for all non-gap letters

		for i in range(al_length):
			targetletter = targetseq[i]
			if targetletter != "-": # meaning anything except gaps
				targetcount += 1
				alignment_column = alignment[:,i] # all letters per site
				nogap_alignment_column = alignment_column.replace("-","").replace("X","") # excluding gaps
				aa_counter = Counter( nogap_alignment_column )

				# calculate float of no-gaps/all characters
				if len(nogap_alignment_column) * 1.0 / len(alignment_column) < gapcutoff:
					# has at least half gaps, for calculate for full length and force less than 50% identity
					conservation = 100.0 * aa_counter[targetletter] / len(alignment_column)
				else: # meaning has more characters than gaps, so calculate normal conservation
					conservation = 100.0 * aa_counter[targetletter] / len(nogap_alignment_column)

				# if reporting raw percentage
				if store_identity:
					index_to_identity[targetcount] = conservation
				else: # otherwise use ranked score from 1 to 9
					if conservation == 100.0:
						rankedscore = 9
					elif conservation >= 98.0:
						rankedscore = 8
					elif conservation >= 95.0:
						rankedscore = 7
					elif conservation >= 90.0:
						rankedscore = 6
					elif conservation >= 80.0:
						rankedscore = 5
					elif conservation >= 70.0:
						rankedscore = 4
					elif conservation >= 60.0:
						rankedscore = 3
					elif conservation >= 50.0:
						rankedscore = 2
					else:
						rankedscore = 1
					index_to_identity[targetcount] = rankedscore
		print >> sys.stderr, "# Calculated conservation for {} sites".format( len(index_to_identity) )
		identindex_dict[target_seqid] = index_to_identity
	return identindex_dict

def rewrite_pdb(pdbfile, seqidlist, scoredict, wayout, forcerecode):
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
		# get relevant chains that match the sequence, in case of hetero multimers
		if record=="DBREF":
			defaultchain = False
			proteinid = line[42:56].strip()
			for seqid in seqidlist:
				if seqid.find(proteinid)>-1:
					chaintarget = line[12]
					print >> sys.stderr, "### keeping chain {} for sequence {}".format( chaintarget, proteinid )
					keepchains[chaintarget] = proteinid
		# DBREF lines should come before ATOM lines, so for all other lines, check for ATOM or not
		if record=="ATOM": # skip all other records
			chain = line[21]
			residue = int( line[22:26] )
			if defaultchain or forcerecode or chain in keepchains: # default chain means take all, or use chain A
				atomcounter += 1
				if defaultchain or forcerecode: # assume only one seqid
					residuescore = scoredict[seqidlist[0]].get(residue,0.00)
				else:
					residuescore = scoredict[keepchains[chain]].get(residue,0.00)
				if residuescore:
					residuecounter[residue] = True
			else: # meaning in another chain, so color as insufficient
				residuescore = 0
			newline = "{}{:6.2f}{}".format( line[:60], residuescore, line[66:].rstrip() )
			print >> wayout, newline
		else: # this will also print DBREF lines
			print >> wayout, line.strip()
	if atomcounter:
		print >> sys.stderr, "# Recoded values for {} atoms in {} residues".format(atomcounter, len(residuecounter) )
	else:
		print >> sys.stderr, "# NO CHAINS FOUND MATCHING SEQ ID {}, CHECK NAME {}".format( seqid, proteinid )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", nargs="*", help="multiple sequence alignment", required=True)
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-g","--gap-cutoff", default=0.5, type=float, help="minimum fraction of non-gap characters per site, else is called unconserved [0.5]")
	parser.add_argument("-i","--identity", action="store_true", help="report percent identity instead of score")
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	parser.add_argument("-s","--sequence", nargs="*", help="sequence ID for PDB", required=True)
	parser.add_argument("--force-recode", action="store_true", help="force recoding regardless of chain")
	args = parser.parse_args(argv)

	if len(args.alignment) != len(args.sequence):
		print >> sys.stderr, "ERROR: {} ALIGNMENTS FOR {} SEQUENCES, MUST BE EQUAL, CHECK -a AND -s".format(len(args.alignment), len(args.sequence)), time.asctime()

	if len(args.sequence) > len(set(args.sequence)):
		print >> sys.stderr, "ERROR: NON UNIQUE NAMES FOR SEQUENCES, CHECK -s"

	conservedict = get_alignment_identity( args.alignment, args.format, args.sequence, args.gap_cutoff, args.identity)
	if conservedict: # indicating that the sequence was found and something was calculated
		rewrite_pdb(args.pdb, args.sequence, conservedict, wayout, args.force_recode)
	else:
		sys.exit("# CANNOT CALCULATE CONSERVATION, EXITING")

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
