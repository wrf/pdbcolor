#!/usr/bin/env python
#
# pdb_heteropecilly.py v1 2017-10-09

'''pdb_heteropecilly.py  last modified 2019-09-25

pdb_heteropecilly.py -a PRP4B_HUMAN.aln -p 4ian.pdb -s PRP4B_HUMAN > 4ian_w_hp.pdb

within a PDB file, fields for atoms are:
Record name      Residue         Position as X Y Z
          Atom serial number                            Occupancy
             Atom    Chain                                   temperatureFactor
                        Residue number
such as:
ATOM      1  N   PRO A  46       0.739  40.031  44.896  1.00 0.842           N

here, the temperatureFactor will be replaced with the heteropecilly decile
as calculated by:

blast_to_align_pairs.py -b simion_taxa/hsapiens_vs_uniprot_blastp.tab -q simion_taxa/Homo_sapiens.fasta.nogaps -s ~/db/human_uniprot.fasta -r simion_taxa/Homo_sapiens.fasta -p heteropecilly-v2/hp_by_site_w_const.tab

residue number must match the alignment, not the position in the model
meaning even if the PDB file starts with residue 20, if the first 19 were
disordered or cleaved, the sequence still must start with residue 1
'''

import sys
import time
import argparse
from collections import Counter
from Bio import AlignIO

def get_alignment_heteropecilly(alignment, alignformat, target_seqid):
	sys.stderr.write("# Reading alignment from {}\n".format( alignment ) )
	alignment = AlignIO.read( alignment, alignformat )

	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)

	sys.stderr.write("# Alignment contains {} taxa for {} sites, including gaps\n".format( num_taxa, al_length ) )

	targetseq = None
	hpscores = None
	for seqrec in alignment:
		try:
			shortid = seqrec.id.split("|")[2]
		except IndexError:
			shortid = None
		if seqrec.id==target_seqid or shortid==target_seqid:
			targetseq = seqrec.seq
		elif seqrec.id=="Heteropecilly_score":
			hpscores = seqrec.seq
	if targetseq is None:
		sys.stderr.write("# ERROR: CANNOT FIND SEQUENCE {}, CHECK OPTION -s OR ALIGNMENT\n".format( target_seqid ) )
		return None
	if hpscores is None:
		sys.stderr.write("# ERROR: CANNOT FIND HETEROPECILLY SCORES {}, CHECK ALIGNMENT\n".format( target_seqid ) )
		return None

	index_to_hp = {}
	nongapcount = 0
	targetcount = 0 # keep track of position in target sequence for all non-gap letters
	for i in range(al_length):
		targetletter = targetseq[i]
		hpvalue = hpscores[i]
		if targetletter != "-": # meaning anything except gaps in target, though there should not be any
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
			index_to_hp[targetcount] = rankedscore
	sys.stderr.write("# Found heteropecilly for {} sites\n".format( nongapcount ) )
	return index_to_hp

def rewrite_pdb(pdbfile, seqid, scoredict, wayout):
	sys.stderr.write("# Reading PDB from {}\n".format(pdbfile) )
	atomcounter = 0
	residuecounter = {}
	keepchains = []
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
			if seqid.find(proteinid)>-1:
				chaintarget = line[12]
				sys.stderr.write("### keeping chain {} for sequence {}\n".format( chaintarget, proteinid ) )
				keepchains.append( chaintarget )
		# for all other lines, check for ATOM or not
		if record=="ATOM": # skip all other records
			chain = line[21]
			residue = int( line[22:26] )
			if defaultchain or chain in keepchains:
				atomcounter += 1
				score = scoredict.get(residue,0.00)
				if score:
					residuecounter[residue] = True
			else: # meaning in another chain, so color as insufficient
				score = 99
			newline = "{}{:6.2f}{}".format( line[:60], score, line[66:].rstrip() )
			print >> wayout, newline
		else:
			print >> wayout, line.strip()
	sys.stderr.write("# Recoded values for {} atoms in {} residues\n".format(atomcounter, len(residuecounter) ) )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", help="multiple sequence alignment", required=True)
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-i","--identity", action="store_true", help="report percent identity instead of score")
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	parser.add_argument("-s","--sequence", help="sequence ID for PDB", required=True)
	args = parser.parse_args(argv)

	hpdict = get_alignment_heteropecilly( args.alignment, args.format, args.sequence)
	if hpdict: # indicating that the sequence was found and something was calculated
		rewrite_pdb(args.pdb, args.sequence, hpdict, wayout)
	else:
		sys.exit("# CANNOT CALCULATE CONSERVATION, EXITING")

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
