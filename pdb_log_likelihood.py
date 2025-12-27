#!/usr/bin/env python
#
# pdb_log_likelihood.py v1 2018-02-01
# v1.1 2019-09-25
# v1.2 2025-12-27 generate PyMOL script

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

def get_alignment_values(alignmentlist, alignformat, targetidlist, print_verbose):
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
			elif seqrec.id in ["Likelihood_score", "Likelihood3_score", "Heteropecilly_score"]:
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
		if scoretype=="Likelihood_score" or scoretype=="Likelihood3_score":
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
	if print_verbose:
		sys.stderr.write("# SCORES ARE: T1:{}, T2:{}, T3:{}\n".format(score_counter[1] + score_counter[2], score_counter[4] + score_counter[5], score_counter[7] + score_counter[8]) )
		for k,v in score_counter.items():
			sys.stderr.write("#{}\t{}\n".format(k,v) )
	return scoreindex_dict

##############################

def get_chains_only(defaultchain, seqidlist, pdbfile):
	'''read PDB file and return two dicts, one where key is chain and value is sequence ID, other where key is the chain and value is integer of the DBREF offset'''
	keepchains = {} # dict where key is chain and value is seqid, though value is not used
	refoffsets = {} # key is chain, value is integer offset from DB seq
	sys.stderr.write("# Reading chain from PDB {}\n".format(pdbfile) )
	for line in open(pdbfile,'r'):
		record = line[0:6].strip()
		# get relevant chains that match the sequence, in case of hetero multimers
		if record=="DBREF":
			defaultchain = False
			proteinid = line[42:56].strip()
			for seqid in seqidlist:
				if seqid.find(proteinid)>-1:
					chaintarget = line[12]
					chainstart = int(line[14:18].strip())
					dbstart = int(line[55:60].strip())
					chainoffset = dbstart - chainstart
					sys.stderr.write("### keeping chain {} for sequence {} with offset {}\n".format( chaintarget, proteinid, chainoffset ) )
					keepchains[chaintarget] = proteinid
					refoffsets[chaintarget] = chainoffset
	if defaultchain: # meaning nothing was found, use default and single sequence
		if seqidlist: # all default chains are assumed to use only the first sequence
			keepchains[defaultchain] = seqidlist[0]
		else: # value is not called, but just to indicate that -s is or used or not
			keepchains[defaultchain] = "UNKNOWN"
		refoffsets[defaultchain] = 0
		sys.stderr.write("### using default chain {}\n".format( defaultchain ) )
	return keepchains, refoffsets

##############################

def make_output_script(keepchains, refoffsets, scoredict, whitebg, wayout):
	'''from the lnl information, print a script for PyMOL'''

	# each color is RGB triplet, values ranging from 0 to 1, where black is 0 and white is 1
	# colors should be:
	# for values [ -1 gray
	#            0 clay   1 pink   2 red
	#            3 mud    4 teal   5 green
	#            6 lead   7 cobalt 8 blue
	#            9 orange ]

	colors = [ [0.39, 0.39, 0.39] ,
             [0.49,0.36,0.40] , [0.76,0.26,0.56] , [0.77,0.20,0.33] ,
             [0.42,0.49,0.42] , [0.26,0.77,0.64] , [0.19,0.74,0.34] ,
             [0.39,0.37,0.55] , [0.35,0.26,0.76] , [0.24,0.45,0.87] ,
             [1.00,0.60,0.12] ]

	if whitebg: # replaces colors -1, 0, 3, 6, 9
		colors = [ [0.99,0.87,0.37] ,
             [0.83,0.74,0.76] , [0.76,0.26,0.56] , [0.77,0.20,0.33] ,
             [0.73,0.82,0.73] , [0.26,0.77,0.64] , [0.19,0.74,0.34] ,
             [0.74,0.73,0.86] , [0.35,0.26,0.76] , [0.24,0.45,0.87] ,
             [0.84,0.46,0.00] ]

	groupnames = [ "gaps", 
                  "t1-weak", "t1-strong", "t1-max",
                  "t2-weak", "t2-strong", "t2-max",
                  "t3-weak", "t3-strong", "t3-max",
                   "const" ]
	group_values = [-1,  0,1,2,  3,4,5,  6,7,8,  9] # encoded earlier

	# begin printing commands for PyMOL script
	wayout.write("hide everything\n")
	if whitebg:
		wayout.write("bg white\n")
	else:
		wayout.write("bg black\n")
	chainstring = " or ".join( ["chain {}".format(chain) for chain in keepchains.keys()] )
	wayout.write("show cartoon, ({})\n".format(chainstring) )

	for i,rgb in enumerate(colors):
		colorname = "{}-col".format( groupnames[i] )
		wayout.write("set_color {}, [{}]\n".format( colorname, ",".join(map(str,rgb)) ) )

	# make commands for each chain
	for chain,seqid in keepchains.items(): # keys are chain letters, values are seq IDs
		chainoffset = refoffsets.get(chain, 0)
		scoregroups = defaultdict(list) # key is group_values group, value is list of residues
		# for each residue, assign to a bin
		for residue in scoredict[seqid].keys():
			residuescore = scoredict[seqid].get(residue,None)
			scoregroups[residuescore].append(residue - chainoffset)
		# for each bin, make a command to color all residues of that bin
		for i,value in enumerate(group_values):
			binname = "{}_{}".format( groupnames[i], chain )
			resilist = list(map(str,scoregroups.get(value,[])))
			if resilist: # do not print empty groups
				binresidues = ",".join(resilist)
				wayout.write("select {}, (chain {} & resi {})\n".format( binname, chain, binresidues ) )
				wayout.write("color {}-col, {}\n".format( groupnames[i], groupnames[i], binname ) )
	# no return

##############################

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

##############################

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", nargs="*", help="multiple sequence alignment", required=True)
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-p","--pdb", help="PDB format file, if none given, then print PyMOL script as output", required=True )
	parser.add_argument("-s","--sequence", nargs="*", help="sequence ID for PDB", required=True )
	parser.add_argument("-w","--write-script", action="store_true", help="write to script instead of recoding PDB file")
	parser.add_argument("-v","--verbose", action="store_true", help="extra output")
	parser.add_argument("--white-bg", action="store_true", help="change colors for white background, i.e. for paper")
	parser.add_argument("--default-chain", default="A", help="default letter of chain when writing to script [A], if DBREF for the sequence cannot be found in PDB")

	parser.add_argument("--color-gaps", action="store_true", help="use dark colors to uniquely color gaps of up to 6 chains")
	parser.add_argument("--force-recode", action="store_true", help="force recoding regardless of chain")
	parser.add_argument("--heteroatoms", action="store_true", help="color heteroatoms, in recoding mode")
	args = parser.parse_args(argv)

	if len(args.alignment) != len(args.sequence):
		sys.stderr.write("ERROR: {} ALIGNMENTS FOR {} SEQUENCES, MUST BE EQUAL, CHECK -a AND -s\n".format(len(args.alignment), len(args.sequence)) )

	scoredict = get_alignment_values( args.alignment, args.format, args.sequence, args.verbose )
	if scoredict: # indicating that the sequence was found and something was calculated
		if args.write_script: # write script to stdout
			refchains, refoffsets = get_chains_only(args.default_chain, args.sequence, args.pdb)
			make_output_script(refchains, refoffsets, scoredict, args.white_bg, wayout)
		else: # assume that user wants to recode the PDB file
			rewrite_pdb(args.pdb, args.sequence, scoredict, wayout, args.force_recode, args.color_gaps, args.heteroatoms)

	else:
		sys.exit("# CANNOT CALCULATE LIKELIHOODS, EXITING")

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
