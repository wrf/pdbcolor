#!/usr/bin/env python
#
# sca_to_pymol_script.py v1

'''sca_to_pymol_script.py v1 last modified 2018-06-20

sca_to_pymol_script.py -a aequorins.aln -e vec_by_site.tab -p 1EJ3.pdb -s Aequorin1

    script is automatically named with .sectors.py
    thus 1EJ3.pdb becomes 1EJ3.pdb.sectors.py

    run in the PyMOL console as:
@1EJ3.pdb.sectors.py
'''

import sys
import time
import os
import argparse
from collections import defaultdict
from Bio import AlignIO

single_letter_conv = { "ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F",
                       "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", "LEU":"L",
                       "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R",
                       "SER":"S", "THR":"T", "VAL":"V", "TRP":"W", "TYR":"Y" }

def read_sca_eigenvalues(eigenfile):
	'''read eigenvalues and return dict where key is sector, value is dict of residues and values'''
	sectors = defaultdict( lambda: defaultdict(float) )
	residuecounter = 0
	sys.stderr.write("# reading eigenvalues from {}".format(eigenfile) + time.asctime() + os.linesep )
	for line in open(eigenfile, 'r'):
		line = line.strip()
		if line and line[1]!="S":
			lsplits = line.replace('"','').split("\t") # remove quote marks
			position = int(lsplits[0])
			e2val = float(lsplits[2])
			e4val = float(lsplits[4])
			if e2val >= 0.05: # blue sector, 1
				residuecounter += 1
				sectors[1][position] = e2val
			elif e2val <= -0.05: # red sector, 2
				residuecounter += 1
				sectors[2][position] = abs(e2val)
			elif e4val >= 0.05: # green sector, 3
				residuecounter += 1
				sectors[3][position] = e4val
	sys.stderr.write("# found {} sectors for {} residues\n".format(len(sectors), residuecounter) )
	return sectors

def calculate_rgb(value, sector):
	'''return list of three RGB values as floats'''
	if value > 0.4: # set max to 0.4
		value = 0.4
	firstcorrector = 0.8 - (value * 2.0)
	secondcorrector = 0.8 - (value * 1.7)
	if sector==1: # blue
		RGB = [firstcorrector, secondcorrector, 0.8]
	elif sector==2: # red
		RGB = [0.8, secondcorrector ,firstcorrector]
	else: # sector==3 # green
		RGB = [secondcorrector, 0.8, firstcorrector]
	return RGB

def get_site_offset(target_seqid, alignfile, alignformat, pdbfile):
	'''read alignment and PDB, get spacing for target sequence, return dict where key is alignment position and value is residue number in structure'''
	sys.stderr.write("# Reading alignment from {}\n".format( alignfile ) )
	alignment = AlignIO.read( alignfile, alignformat )
	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)
	sys.stderr.write("# Alignment contains {} taxa for {} sites, including gaps\n".format( num_taxa, al_length ) )
	targetseq = None
	for seqrec in alignment:
		if seqrec.id==target_seqid:
			targetseq = seqrec.seq
	if targetseq is None:
		sys.stderr.write("# ERROR: CANNOT FIND SEQUENCE {}, CHECK OPTION -s OR ALIGNMENT\n".format( target_seqid ) )
		return None
	target_position = 0
	aapos_to_align = {} # key is prot position, value is tuple of position in alignment and amino acid
	aapos_to_res = {} # key is prot position, value is amino acid at that position
	sys.stderr.write("# Determining alignment positions for {}\n".format(target_seqid) + time.asctime() + os.linesep )
	for i in range(al_length):
		targetletter = targetseq[i] # seq index is 0 to max
		if targetletter != "-": # meaning anything except gaps
			target_position += 1
			aapos_to_align[target_position] = i+1 # i starts at 0, must offset to 1 for residue index
			aapos_to_res[target_position] = targetletter
	total_residues = len(aapos_to_res)

	residue_offset = 0 # start at 0, meaning first residue in alignment is first in protein
	advance_offset = True # one time flag to offset based on first residue in structure
	alignpos_to_res = {} # key is alignment position, value is residue number in structure
	found_residues = {}
	#ATOM      1  N   ASP A   3     -32.800  10.288   5.960  1.00 30.00           N
	sys.stderr.write("# Reading protein structure from {}\n".format( pdbfile ) )
	for line in open(pdbfile, 'r'):
		record = line[0:6].strip()
		if record=="ATOM": # skip all other records
			chain = line[21]
			residue = single_letter_conv.get(line[17:20], None)
			if residue is None:
				continue
			pdb_residue_num = int( line[22:26] )
			if pdb_residue_num in found_residues: # do not process if this residue was found for another ATOM
				continue
			if advance_offset: # if first residue in structure is ahead of alignment
				sys.stderr.write("# First residue {} in structure is {}, offsetting\n".format(residue, pdb_residue_num) )
				advance_offset = False # stop checking after first time
				residue_offset = (pdb_residue_num-1)
			found_residues[pdb_residue_num] = True
			align_res_num = pdb_residue_num + residue_offset
			#sys.stderr.write(residue,  pdb_residue_num, residue_offset, align_res_num, aapos_to_res.get(align_res_num,"X"), aapos_to_align[align_res_num])
			if residue==aapos_to_res.get(align_res_num,"X"):
				alignpos_to_res[aapos_to_align[align_res_num]] = pdb_residue_num
			else: # residue is offset, so increment
				residue_offset -= 1
				#sys.stderr.write("offset = {}".format(residue_offset))
	return alignpos_to_res

def make_pymol_script(pdbfile, protsectors, aligntores):
	'''generate a pymol script with coloring instructions'''
	script_file = "{}.sectors.py".format(pdbfile)
	sys.stderr.write("# Generating python script {}\n".format( script_file ) )
	with open(script_file,'w') as sf:
		sf.write("bg white\n")
		sf.write("hide all\n")
		sf.write("show cartoon\n")
		sf.write("set_color colordefault, [0.8,0.8,0.8]\n")
		sf.write("color colordefault, all\n")
		for sector,residict in protsectors.items():
			sectorname = "sector_{}".format(sector)
			resilist = [str(aligntores.get(k,None)) for k in sorted(residict.keys()) if not None]
			sectorresidues = ",".join(resilist)
			sf.write("select {}, (resi {})\n".format( sectorname, sectorresidues ) )
			sf.write("show spheres, {}\n".format(sectorname) )
			for residue in sorted(residict.keys()):
				adjresidue = aligntores.get(residue,None)
				if adjresidue is not None:
					sf.write("set_color colres{}, [{},{},{}]\n".format(adjresidue, *calculate_rgb(residict[residue], sector)) )
					sf.write("color colres{}, (resi {})\n".format( adjresidue, adjresidue ) )
	sys.stderr.write("# Finished writing python script {}".format( script_file ) + time.asctime() + os.linesep)

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", help="multiple sequence alignment", required=True)
	parser.add_argument("-e","--eigenvalues", help="tabular file of eigenvalues by site from R script")
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	parser.add_argument("-s","--sequence", help="sequence ID for PDB", required=True)
	args = parser.parse_args(argv)

	siteoffsetdict = get_site_offset(args.sequence, args.alignment, args.format, args.pdb)
	sectordict = read_sca_eigenvalues(args.eigenvalues)
	make_pymol_script(args.pdb, sectordict, siteoffsetdict)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
