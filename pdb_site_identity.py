#!/usr/bin/env python
#
# pdb_site_identity.py v1 2017-07-25

'''pdb_site_identity.py  last modified 2022-12-28

    two usage options: make a PyMol script (-w), or rewrite the PDB file

  proteins or ligands not in the alignment are coded as -1, and colored "null"

  if recoding the PDB file:

pdb_site_identity.py -a mox_all.aln -s DOPO_HUMAN -p 4zel.pdb > 4zel_w_scores.pdb

within a PDB file, fields for atoms are:
Record name      Residue         Position as X Y Z
          Atom serial number                            Occupancy
             Atom    Chain                                   temperatureFactor
                        Residue number
such as:
ATOM      1  N   PRO A  46       0.739  40.031  44.896  1.00 0.842           N

    the temperatureFactor will be replaced with the identity %,
    where bins are:
  0, 50, 60, 70, 80, 90, 95, 98, 100%

residue number must match the alignment, not the position in the model
meaning even if the PDB file starts with residue 20, if the first 19 were
disordered or cleaved, the sequence still must start with residue 1

meaning it is important that the alignment is not trimmed or truncated

this is also important if 6-His tags or similar are in the structure, the
offset can be changed using the option -I / --ignore-offset

most protein models (like AlphaFold) will start with residue 1
'''

import sys
import os
import time
import argparse
from collections import Counter,defaultdict
from numpy import log
from Bio import AlignIO

def get_conservation(alignmentlist, alignformat, targetidlist):
	'''read alignment and calculate conservation scores for each position, and return a dict where key is seqid, value is a dict of positions and scores

    positional conservation scores defined by:
    D (divergence) = f(i,a) * ln( f(i,a)/q(a) ) + (1 - f(i,a)) * ln ( (1-f(i,a))/(1-q(a)) )
    f(i,a) is frequency of amino acid 'a' at position 'i'
    q(a) is background frequency of amino acid 'a' in all proteins in the alignment'''
	consindex_dict = {} # key is target seqid, value is dict of conservation scores
	for alignment, target_seqid in zip(alignmentlist,targetidlist):
		sys.stderr.write("# Reading alignment from {}\n".format( alignment ) )
		alignment = AlignIO.read( alignment, alignformat )

		al_length = alignment.get_alignment_length()
		num_taxa = len(alignment)

		sys.stderr.write("# Alignment contains {} taxa for {} sites, including gaps\n".format( num_taxa, al_length ) )

		targetseq = None
		for seqrec in alignment:
			if seqrec.id==target_seqid:
				targetseq = seqrec.seq
		if targetseq is None:
			sys.stderr.write("ERROR: CANNOT FIND SEQUENCE {}, CHECK OPTION -s OR ALIGNMENT\n".format( target_seqid ) )
			return None

		index_to_cons = {}
		targetcount = 0 # keep track of position in target sequence for all non-gap letters

		base_aa_counts = Counter()
		sys.stderr.write("# Calculating global amino acid frequencies  {}\n".format( time.asctime() ) )
		for i in range(al_length): # only iterate over columns that are not gaps in target seq
			targetletter = targetseq[i]
			if targetletter != "-": # meaning anything except gaps
				alignment_column = alignment[:,i] # all letters per site
				base_aa_counts.update(alignment_column)
		total_base_aas = sum(base_aa_counts.values()) - base_aa_counts.get("-",0)
		base_frequencies = dict([ (AA, base_aa_counts[AA]*1.0/total_base_aas) for AA in "ACDEFGHIKLMNPQRSTVWY-"])
		#sys.stderr.write(base_frequencies)

		sys.stderr.write("# Calculating sitewise conservation  {}\n".format( time.asctime() ) )
		for i in range(al_length):
			targetletter = targetseq[i]
			if targetletter != "-": # meaning anything except gaps
				targetcount += 1
				alignment_column = alignment[:,i] # all letters per site
				nogap_alignment_column = alignment_column.replace("-","").replace("X","") # excluding gaps
				aa_counter = Counter( nogap_alignment_column )

				# calculate float of no-gaps/all characters
				target_freq = aa_counter[targetletter] * 1.0 / len(nogap_alignment_column)
				if target_freq == 1: # otherwise will divide by zero at 1-f(i,a), so simplifies
					conservation = log(1/base_frequencies[targetletter])
				else: # calculate by formula
					conservation = target_freq * log(target_freq/base_frequencies[targetletter]) + ( (1-target_freq) * log( (1-target_freq)/(1-base_frequencies[targetletter]) ) )
				# if reporting raw percentage
				index_to_cons[targetcount] = conservation
		sys.stderr.write("# Calculated conservation for {} sites\n".format( len(index_to_cons) ) )
		consindex_dict[target_seqid] = index_to_cons
	return consindex_dict

def get_alignment_identity(alignmentlist, alignformat, targetidlist, gapcutoff):
	'''read alignment from file, return a dict where key is seqid, value is a dict where keys are integers of position of the target seq in the alignment (starting from 1), and value is float of identity'''
	identindex_dict = {} # key is target seqid, value is dict of identities
	for alignment, target_seqid in zip(alignmentlist,targetidlist):
		sys.stderr.write("# Reading alignment from {}\n".format( alignment ) )
		alignment = AlignIO.read( alignment, alignformat )

		al_length = alignment.get_alignment_length()
		num_taxa = len(alignment)

		sys.stderr.write("# Alignment contains {} taxa for {} sites, including gaps\n".format( num_taxa, al_length ) )

		targetseq = None
		for seqrec in alignment:
			if seqrec.id==target_seqid:
				targetseq = seqrec.seq
		if targetseq is None:
			sys.stderr.write("ERROR: CANNOT FIND SEQUENCE {} IN ALIGNMENT, CHECK OPTION -s OR ALIGNMENT\n".format( target_seqid ) )
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
				index_to_identity[targetcount] = conservation
		sys.stderr.write("# Calculated sitewise identity for {} sites\n".format( len(index_to_identity) ) )
		identindex_dict[target_seqid] = index_to_identity
	return identindex_dict


def parse_dbref_line(dbref_line, seqidlist, forcerecode, ignoreoffset):
	keepchains = {} # dict where key is chain and value is seqid
	refoffsets = {} # key is chain, value is integer offset from DB seq
	proteinid = dbref_line[42:56].strip()
	chaintarget = dbref_line[12]
	sys.stderr.write("# found DBREF info for seq {} as chain {}\n".format( proteinid , chaintarget ) )
	for seqid in seqidlist:
		if seqid.find(proteinid)>-1 or forcerecode:
			chainstart = int(dbref_line[14:18].strip())
			dbstart = int(dbref_line[55:60].strip())
			chainoffset = dbstart - chainstart
			if forcerecode:
				sys.stderr.write("### forcing recode on chain {} for seq {} with seq {}, starting at {} with offset {}\n".format( chaintarget, proteinid, seqidlist[0], chainstart, chainoffset ) )
				keepchains[chaintarget] = seqidlist[0]
			else:
				sys.stderr.write("# keeping chain {} for {}, starting at {} with offset {}\n".format( chaintarget, proteinid, chainstart, chainoffset ) )
				keepchains[chaintarget] = proteinid
			if ignoreoffset is not None:
				chainoffset = ignoreoffset
				sys.stderr.write("### forcing offset to {}\n".format(chainoffset) )
			refoffsets[chaintarget] = [chainstart, chainoffset]
	else: # meaning no seq found
		pass
		#sys.stderr.write("### NOTE: no seq in aligment matches seq {}\n".format( proteinid ) )
	return keepchains, refoffsets


def get_chains_only(defaultchain, seqidlist, pdbfile, forcerecode, ignoreoffset):
	'''read PDB file and return two dicts, one where key is chain and value is sequence ID, other where key is the chain and value is integer of the DBREF offset'''
	keepchains = {} # dict where key is chain and value is seqid
	refoffsets = {} # key is chain, value is integer offset from DB seq
	sys.stderr.write("# Reading chain from PDB {}\n".format(pdbfile) )
	for line in open(pdbfile,'r'):
		record = line[0:6].strip()
		# get relevant chains that match the sequence, in case of hetero multimers
		if record=="DBREF":
			defaultchain = False # change to False, indicating that chains were found
			line_kc_d, line_ro_d = parse_dbref_line(line, seqidlist, forcerecode, ignoreoffset)
			keepchains.update(line_kc_d)
			refoffsets.update(line_ro_d)
	if defaultchain: # meaning nothing for DBREF was found, use default and single sequence
		sys.stderr.write("# WARNING: NO DBREF TAGS FOUND {} , check PDB file or force recoding with --force-recode \n".format( " ".join(seqidlist) ) )
		keepchains[defaultchain] = seqidlist[0]
		refoffsets[defaultchain] = [1,0]
	return keepchains, refoffsets


def rewrite_pdb(pdbfile, seqidlist, scoredict, wayout, forcerecode, ignoreoffset):
	sys.stderr.write("# Reading PDB from {}\n".format(pdbfile) )
	atomcounter = 0
	residuecounter = {}
	keepchains = {} # dict where key is chain and value is seqid
	refoffsets = {} # key is chain, value is integer offset from DB seq
	defaultchain = True # flag for whether DBREF occurs at all
	for line in open(pdbfile,'r'):
		# records include:
		# HEADER TITLE COMPND SOURCE AUTHOR REVDAT JRNL REMARK
		# DBREF SEQRES HET HETNAM FORMUL HELIX SHEET SSBOND LINK CISPEP SITE ATOM CONECT

		record = line[0:6].strip()
		# get relevant chains that match the sequence, in case of hetero multimers
		if record=="DBREF":
			defaultchain = False
			line_kc_d, line_ro_d = parse_dbref_line(line, seqidlist, forcerecode, ignoreoffset)
			keepchains.update(line_kc_d)
			refoffsets.update(line_ro_d)
		# DBREF lines should come before ATOM lines, so for all other lines, check for ATOM or not
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

		if record=="ATOM": # skip all other records
			chain = line[21]
			residue = int( line[22:26] )
			chainstart, chainoffset = refoffsets.get(chain, [1,0] )
			if residue < chainstart:
				if residue < 1: # chain includes expression vector, 6-HIS tag or similar
					sys.stderr.write("# SKIPPING NEGATIVE RESIDUE {} (EXP VECTOR)\n".format(residue) )
				continue
			if defaultchain or forcerecode or chain in keepchains: # default chain means take all, or use chain A
				atomcounter += 1
				if defaultchain or forcerecode: # assume only one seqid
					residuescore = scoredict[seqidlist[0]].get(residue+chainoffset,0.00)
				else:
					residuescore = scoredict[keepchains[chain]].get(residue+chainoffset,0.00)
				if residuescore:
					residuecounter[residue] = True
			else: # meaning in another chain, so color as insufficient
				residuescore = -1
			newline = "{}{:6.2f}{}\n".format( line[:60], residuescore, line[66:].rstrip() )
			wayout.write(newline)
		else: # this will also print DBREF lines
			wayout.write(line)
	if atomcounter:
		sys.stderr.write("# Recoded values for {} atoms in {} residues\n".format(atomcounter, len(residuecounter) ) )
	else:
		sys.stderr.write("# ERROR: NO CHAINS FOUND MATCHING SEQ ID {}\n".format( seqid ) )


def make_output_script(scoredict, keepchains, refoffsets, colorscript, forcerecode, basecolor="red"):
	'''from the identity calculations, print a script for PyMOL'''
	colors = {} # key is colorscheme name, value is list of colors
	colors["red"] = [ [0.63,0.63,0.63] , [0.73,0.55,0.55] , [0.75,0.47,0.47], 
				   [0.77,0.38,0.38] , [0.79,0.29,0.29] , [0.82,0.21,0.21], 
				   [0.84,0.13,0.13]    , [0.88,0,0]     , [1,0,0.55] ]
	colors["yellow"] = [ [0.63,0.63,0.63] , [0.66,0.67,0.51] , [0.68,0.71,0.43], 
				   [0.70,0.73,0.36] , [0.72,0.76,0.28] , [0.74,0.80,0.20], 
				   [0.76,0.83,0.12] , [0.79,0.87,0.00] , [1,0.8,0.16] ]
	colors["blue"] = [ [0.63,0.63,0.63] , [0.50,0.58,0.68] , [0.42,0.55,0.71], 
				   [0.35,0.52,0.73] , [0.28,0.49,0.76] , [0.20,0.46,0.80], 
				   [0.12,0.43,0.83] , [0.00,0.38,0.87] , [0.6,0,1] ]
	colors["green"] = [ [0.63,0.63,0.63] , [0.50,0.68,0.56] , [0.42,0.71,0.53], 
				   [0.35,0.74,0.49] , [0.26,0.77,0.44] , [0.19,0.80,0.41], 
				   [0.12,0.83,0.37] , [0.01,0.87,0.31] , [0,1,0.83] ]
	insuf_color = [0.75, 0.75, 0.58]

	binvalues = [0.0, 50.0 ,60.0 ,70.0 ,80.0 ,90.0 ,95.0 ,98.0 ,100, 101]

	if len(keepchains) < 1: # meaning no chains kept
		sys.stderr.write("# ERROR: no chains found, cannot generate script {} , check PDB or force recode\n".format(colorscript) )
		return None

	sys.stderr.write("# Generating PyMOL script {}\n".format(colorscript) )
	# begin printing commands for PyMOL script
	with open(colorscript, 'w') as cs:
		cs.write("hide everything\n")
		cs.write("show cartoon\n")
		cs.write("set_color colordefault, [{}]\n".format( ",".join(map(str,insuf_color)) ) )
		cs.write("color colordefault, all\n" )
		# make commands for each color
		for color, colorlist in colors.items():
			for i,rgb in enumerate(colorlist):
				colorname = "{}{}".format( color, int(binvalues[i]) )
				cs.write("set_color {}, [{}]\n".format( colorname, ",".join(map(str,rgb)) ) )

		# make commands for each chain
		for chain in keepchains.keys(): # keys are chain letters, values are seq IDs
			chainstart, chainoffset = refoffsets.get(chain, [1,0] )
			pctgroups = defaultdict(list) # key is percent group, value is list of residues
			# for each residue, assign to a bin
			for residue in scoredict[keepchains[chain]].keys():
				if residue < chainstart: # chain begins partway through protein
					if residue < 1: # chain includes expression vector, 6-HIS tag or similar
						sys.stderr.write("# SKIPPING NEGATIVE RESIDUE {} (EXP VECTOR)\n".format(residue) )
					continue
				residuescore = scoredict[keepchains[chain]].get(residue,0.00)
				# determine which bin the residue belongs in
				for i,value in enumerate(binvalues[:-1]): # iterate through bins
					upper_bin_limit = binvalues[i+1] # upper limit would be next bin up
					if residuescore < upper_bin_limit: # check if residue can be in that bin, otherwise try next bin up
						adjusted_residue = residue - chainoffset
						if adjusted_residue > 0: # skip negative residues, meaning offset is big and structure starts after actual gene
							pctgroups[value].append(adjusted_residue)
							break # once bin is found, break loop
				# should not need an else
			# assign whole chain to lowest color, then build up
			cs.write("color {}0, chain {}\n".format( basecolor, chain ) )
			# for each bin, make a command to color all residues of that bin
			for i,value in enumerate(binvalues[:-1]):
				if i==0: # long lists apparently crash the program, so skip
					continue
				binname = "{}pct_grp_{}_{}".format( int(value), i+1, chain )
				resilist = list(map(str,pctgroups[value]))
				if len(resilist)==0: # empty list
					continue
				binresidues = ",".join(resilist)
				cs.write("select {}, (chain {} & resi {})\n".format( binname, chain, binresidues ) )
				cs.write("color {}{}, {}\n".format( basecolor, int(value), binname ) )
	sys.stderr.write("# Run as:\n@{}\n".format( os.path.abspath(colorscript) ) )
	return colorscript

def print_stats(identitydict):
	'''based on the identity scores, print a short table indicating overall conservation'''
	for seqid,posscores in identitydict.items():
		num_sites = len(posscores.values())
		averageid = sum(posscores.values()) * 1.0 / num_sites
		lowestid = min(posscores.values())
		highestid = max(posscores.values())
		num_identities = list(posscores.values()).count(highestid)
		sys.stderr.write("{}\tlength:{}\tmean:{:.2f}\tmin:{:.1f}\tmax:{:.1f}\tN_max:{} ({:.1f}%)\n".format( seqid, num_sites, averageid, lowestid, highestid, num_identities, num_identities*100.0/num_sites) )
	# no return

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", nargs="*", help="multiple sequence alignment : required", required=True)
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-g","--gap-cutoff", metavar="FLOAT from 0.0-1.0", default=0.5, type=float, help="minimum fraction of non-gap characters per site, else is called unconserved [default is 0.5]")
	parser.add_argument("-p","--pdb", metavar="FILENAME", help="PDB format file : required", required=True)
	parser.add_argument("-s","--sequence", nargs="*", help="sequence ID for PDB, protein seq ID that matches the DBREF field in the PDB file : required", required=True)
	parser.add_argument("-w","--write-script", metavar="OUTPUT-SCRIPT-NAME", help="write to script instead of recoding PDB file")
	parser.add_argument("--base-color", default="red", help="color gradient when writing to script, default is red, options are [red,yellow,green,blue]")
	parser.add_argument("--default-chain", default="A", help="default letter of chain when writing to script [A], if DBREF for the sequence cannot be found in PDB")
	parser.add_argument("--ct-conservation", action="store_true", help="calculate sitewise positional conservation (CT model)")
	parser.add_argument("--force-recode", action="store_true", help="force recoding regardless of chain, if DBREF is missing or the protein name does not match the alignment")
	parser.add_argument("-I", "--ignore-offset", nargs='?', type=int, const=0, help="do not read offset from PDB DBREF, use 0, or specify integer")
	parser.add_argument("--stats", action="store_true", help="print basic stats")
	args = parser.parse_args(argv)

	if len(args.alignment) != len(args.sequence):
		sys.stderr.write("ERROR: {} ALIGNMENTS FOR {} SEQUENCES, MUST BE EQUAL, CHECK -a AND -s  {}\n".format(len(args.alignment), len(args.sequence), time.asctime() ) )

	if len(args.sequence) > len(set(args.sequence)):
		sys.stderr.write("ERROR: NON UNIQUE NAMES FOR SEQUENCES, CHECK -s\n")

	# calculate identity or conservation
	if args.ct_conservation:
		conservedict = get_conservation( args.alignment, args.format, args.sequence)
	else:
		conservedict = get_alignment_identity( args.alignment, args.format, args.sequence, args.gap_cutoff)

	# print stats
	if args.stats:
		print_stats(conservedict)

	# rewrite PDB file or make PyMOL script
	if conservedict: # indicating that the sequence was found and something was calculated
		if args.write_script: # write output PyMOL script with color commands
			refchains, refoffsets = get_chains_only(args.default_chain, args.sequence, args.pdb, args.force_recode, args.ignore_offset)
			output_script_command = make_output_script(conservedict, refchains, refoffsets, args.write_script, args.force_recode, args.base_color)
		else: # recode PDB beta-factors
			rewrite_pdb(args.pdb, args.sequence, conservedict, wayout, args.force_recode, args.ignore_offset)
	else:
		sys.exit("# CANNOT CALCULATE CONSERVATION, EXITING")

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
