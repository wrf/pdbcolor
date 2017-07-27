#!/usr/bin/env python
#
# pdb_site_identity.py v1 2017-07-25

'''pdb_site_identity.py  last modified 2017-07-25

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

'''

import sys
import time
import argparse
from collections import Counter
from Bio import AlignIO

def get_alignment_conservation(alignment, alignformat, target_seqid, gapcutoff, store_identity=False):
	print >> sys.stderr, "# Reading alignment from {}".format( alignment )
	alignment = AlignIO.read( alignment, alignformat )

	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)

	print >> sys.stderr, "# Alignment contains {} taxa for {} sites".format( num_taxa, al_length )

	targetseq = None
	for seqrec in alignment:
		if seqrec.id==target_seqid:
			targetseq = seqrec.seq
	if targetseq is None:
		print >> sys.stderr, "# ERROR: CANNOT FIND SEQUENCE {}".format( target_seqid )

	index_to_identity = {}
	targetcount = 0 # keep track of position in target sequence for all non-gap letters

	for i in range(al_length):
		targetletter = targetseq[i]
		if targetletter != "-": # meaning anything except gaps
			targetcount += 1
			alignment_column = alignment[:,i] # all letters per site
			nogap_alignment_column = alignment_column.replace("-","") # excluding gaps
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
	return index_to_identity

def rewrite_pdb(pdbfile, conservedict, wayout):
	print >> sys.stderr, "# Writing new PDB from {}".format(pdbfile)
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
		if record=="ATOM": # skip all other records
			residue = int( line[22:26] )
			conservescore = conservedict.get(residue,0.00)
			newline = "{}{:6.2f}{}".format( line[:60], conservescore, line[66:].rstrip() )
			print >> wayout, newline
		else:
			print >> wayout, line.strip()

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser()
	parser.add_argument("-a","--alignment", help="multiple sequence alignment", required=True)
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-g","--gap-cutoff", default=0.5, type=float, help="minimum fraction of non-gap characters per site, else is called unconserved [0.5]")
	parser.add_argument("-i","--identity", action="store_true", help="report percent identity instead of score")
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	parser.add_argument("-s","--sequence", help="sequence ID for PDB", required=True)
	args = parser.parse_args(argv)

	conservedict = get_alignment_conservation( args.alignment, args.format, args.sequence, args.gap_cutoff, args.identity)
	rewrite_pdb(args.pdb, conservedict, wayout)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
