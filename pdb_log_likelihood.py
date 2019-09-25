#!/usr/bin/env python
#
# pdb_log_likelihood.py v1 2018-02-01

'''pdb_log_likelihood.py  last modified 2019-09-25

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
import argparse
from collections import Counter,defaultdict
from Bio import AlignIO

def get_alignment_values(alignmentlist, alignformat, targetidlist, printw):
	scoreindex_dict = {} # dict of dicts, where key is seqID, value is dict where key is position
	score_counter = defaultdict(int)
	for alignment,target_seqid in zip(alignmentlist,targetidlist):
		sys.stderr.write("# Reading alignment from {}\n".format( alignment ) )
		alignment = AlignIO.read( alignment, alignformat )

		al_length = alignment.get_alignment_length()
		num_taxa = len(alignment)

		sys.stderr.write("# Alignment contains {} sequences for {} sites, including gaps\n".format( num_taxa, al_length ) )

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
			sys.stderr.write("# ERROR: CANNOT FIND SEQUENCE {}, CHECK OPTION -s OR ALIGNMENT\n".format( target_seqid ) )
			return None
		if sitescores is None:
			sys.stderr.write("# ERROR: CANNOT FIND SITEWISE SCORES {}, CHECK ALIGNMENT\n".format( target_seqid ) )
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
					if lnlvalue == "-": # gaps
						rankedscore = -1
					elif lnlvalue == "x" or lnlvalue == "X": # constant sites
						nongapcount += 1
						rankedscore = 9 # is 9 in T3 version
						#rankedscore = 16 # is 16 in T2 version
					else:
						nongapcount += 1
						rankedscore = int(lnlvalue,16) # base 16 string to integer
					index_to_score[targetcount] = rankedscore
					score_counter[rankedscore] += 1
			sys.stderr.write("# Found likelihood for {} sites for {}\n".format( nongapcount, target_seqid ) )
			scoreindex_dict[target_seqid] = index_to_score
		elif scoretype=="Heteropecilly_score":
			for i in range(al_length):
				targetletter = targetseq[i]
				hpvalue = sitescores[i]
				if targetletter != "-":
					targetcount += 1
					if hpvalue == "-": # gap
						rankedscore = -1
					elif hpvalue == "c": # semi-constant site
						nongapcount += 1
						rankedscore = 10
					elif hpvalue == "C": # constant site
						nongapcount += 1
						rankedscore = 11
					else:
						nongapcount += 1
						rankedscore = int(hpvalue)
					index_to_score[targetcount] = rankedscore
					score_counter[rankedscore] += 1
			sys.stderr.write("# Found heteropecilly for {} sites for {}\n".format( nongapcount, target_seqid ) )
			scoreindex_dict[target_seqid] = index_to_score
	if printw:
		sys.stderr.write("# SCORES ARE: T1:{}, T2:{}, T3:{}\n".format(score_counter[1] + score_counter[2], score_counter[4] + score_counter[5], score_counter[7] + score_counter[8]) )
		for k,v in score_counter.items():
			sys.stderr.write("#{}\t{}\n".format(k,v) )
	return scoreindex_dict

def rewrite_pdb(pdbfile, seqidlist, scoredict, wayout, forcerecode, colorgaps, heterocolors):
	sys.stderr.write("# Reading PDB from {}\n".format(pdbfile) )
	atomcounter = 0
	hetatmcounter = 0
	residuecounter = {} # keys are strings of chain + residue
	keepchains = {} # dict where key is chain and value is seqid
	defaultchain = True # flag for whether DBREF occurs at all
	for line in open(pdbfile,'r'):
		# from http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
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
		if record=="DBREF": # find which chains match seq ID
			defaultchain = False
			proteinid = line[42:56].strip()
			for seqid in seqidlist:
				if seqid.find(proteinid)>-1:
					chaintarget = line[12]
					chainstart = int(line[14:18].strip())
					dbstart = int(line[55:60].strip())
					chainoffset = dbstart - chainstart
					sys.stderr.write("### keeping chain {} for sequence {}, starting at {} with offset {}\n".format( chaintarget, proteinid, chainstart, chainoffset ) )
					keepchains[chaintarget] = proteinid
		# for all other lines, check for ATOM or not
		if record=="ATOM": # skip all other records
			chain = line[21]
			residue = int( line[22:26] )
			if defaultchain or forcerecode or chain in keepchains:
				atomcounter += 1
				if defaultchain or forcerecode: # assume only one seqid
					score = scoredict[seqidlist[0]].get(residue,-1.00)
				else:
					score = scoredict[keepchains[chain]].get(residue,-1.00)
				if score:
					chain_residue = "{}/{}".format(chain, residue)
					residuecounter[chain_residue] = True
			else: # meaning in another chain, so color as null
				score = 99
			newline = "{}{:6.2f}{}".format( line[:60], score, line[66:] )
			wayout.write( newline )

		#COLUMNS       DATA  TYPE     FIELD         DEFINITION
		#-----------------------------------------------------------------------
		# 1 - 6        Record name    "HETATM"
		# 7 - 11       Integer        serial        Atom serial number.
		#13 - 16       Atom           name          Atom name.
		#17            Character      altLoc        Alternate location indicator.
		#18 - 20       Residue name   resName       Residue name.
		#22            Character      chainID       Chain identifier.
		#23 - 26       Integer        resSeq        Residue sequence number.
		#27            AChar          iCode         Code for insertion of residues.
		#31 - 38       Real(8.3)      x             Orthogonal coordinates for X.
		#39 - 46       Real(8.3)      y             Orthogonal coordinates for Y.
		#47 - 54       Real(8.3)      z             Orthogonal coordinates for Z.
		#55 - 60       Real(6.2)      occupancy     Occupancy.
		#61 - 66       Real(6.2)      tempFactor    Temperature factor.
		#77 - 78       LString(2)     element       Element symbol; right-justified.
		#79 - 80       LString(2)     charge        Charge on the atom.

		elif record=="HETATM":
			if heterocolors:
				# colors consist of:
				# bluewhite for oxygen in H2O, as [0.85,0.85,1.00]
				# paleblue for metals, as [0.75,0.75,1.0]
				# brightorange for most other ligands, as [1.00,0.70,0.20]
				heteroname = line[17:20].strip()
				element = line[76:78].strip()
				if heteroname=="HOH":
					score = 21
				elif element=="NA" or element=="K":
					score = 22
				elif element=="MG" or element=="CA" or element=="ZN" or element=="FE":
					score = 23
				else:
					score = 24
			else: # color as null
				score = 99
			newline = "{}{:6.2f}{}".format( line[:60], score, line[66:] )
			wayout.write( newline )
		else:
			wayout.write( line )
	sys.stderr.write("# Recoded values for {} atoms in {} residues\n".format(atomcounter, len(residuecounter) ) )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", nargs="*", help="multiple sequence alignment", required=True)
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	parser.add_argument("-s","--sequence", nargs="*", help="sequence ID for PDB", required=True)
	parser.add_argument("-w","--w", action="store_true", help="extra output")
	parser.add_argument("--color-gaps", action="store_true", help="use dark colors to uniquely color gaps of up to 6 chains")
	parser.add_argument("--force-recode", action="store_true", help="force recoding regardless of chain")
	parser.add_argument("--heteroatoms", action="store_true", help="color heteroatoms")
	args = parser.parse_args(argv)

	if len(args.alignment) != len(args.sequence):
		sys.stderr.write("ERROR: {} ALIGNMENTS FOR {} SEQUENCES, MUST BE EQUAL, CHECK -a AND -s\n".format(len(args.alignment), len(args.sequence)) )

	scoredict = get_alignment_values( args.alignment, args.format, args.sequence, args.w)
	if scoredict: # indicating that the sequence was found and something was calculated
		rewrite_pdb(args.pdb, args.sequence, scoredict, wayout, args.force_recode, args.color_gaps, args.heteroatoms)
	else:
		sys.exit("# CANNOT CALCULATE LIKELIHOODS, EXITING")

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
