#!/usr/bin/env python
#
# gff_cds_to_pymol_script.py v1 2018-06-04

'''gff_cds_to_pymol_script.py  last modified 2018-06-04

gff_cds_to_pymol_script.py -g ppyr_00001.gff -i PPYR_00001-PA > ppyr_00001_color_5dv9.pml

run in PyMOL as (meaning do not use the run command):
@ppyr_00001_color_5dv9.pml
'''

import sys
import time
import argparse

def make_output_script(cds_to_residues, out_of_phase_list, chain, wayout):
	'''from the GFF information, print a script for PyMOL'''

	# colors in order of:
	# blue, green, red, orange, purple
	exon_colors = { 1 :[0.65,0.81,0.89], 2 :[0.12,0.47,0.71], 
                    3 :[0.70,0.87,0.54], 4 :[0.20,0.63,0.17], 
                    5 :[0.98,0.60,0.60], 6 :[0.89,0.10,0.11], 
                    7 :[0.99,0.75,0.44], 8 :[1.00,0.50,0.00], 
                    9 :[0.79,0.70,0.84], 10 :[0.42,0.24,0.60] }
	# { 1:"#a6cee3", 2:"#1f78b4", 3:"#b2df8a", 4:"#33a02c", 5:"#fb9a99", 
	# 6:"#e31a1c", 7:"#fdbf6f", 8:"#ff7f00", 9:"#cab2d6", 10:"#6a3d9a" }
	outofphase_color = [0.3,0.3,0.3] # dark grey

	# begin printing commands for PyMOL script
	print >> wayout, "hide everything"
	print >> wayout, "bg white"
	print >> wayout, "show cartoon, (chain {})".format(chain)
	for exonnum, residues in enumerate(cds_to_residues):
		colorlist = [ str(i) for i in exon_colors[exonnum+1] ]
		print >> wayout, "set_color excolor_{}, [{}]".format( exonnum+1, ",".join(colorlist) )
		print >> wayout, "select exon_{}, (resi {}) & (chain {})".format( exonnum+1, ",".join(residues), chain )
		print >> wayout, "color excolor_{}, exon_{}".format( exonnum+1, exonnum+1 ) 
	outofphase_col_str = [ str(i) for i in outofphase_color ]
	print >> wayout, "set_color outphase, [{}]".format( ",".join(outofphase_col_str) )
	print >> wayout, "select out_of_phase, (resi {}) & (chain {})".format( ",".join(out_of_phase_list), chain )
	print >> wayout, "color outphase, out_of_phase"
	# no return

def attributes_to_dict(attributes):
	'''convert GFF attribute string into dictionary of key-value pairs'''
	attrd = {}
	#if attributes.find("ID")==0: # indicates GFF3 format
	# if one of the terms does not have = sign, perhaps Note, then ignore
	attrd = dict([(field.strip().split("=")) for field in attributes.split(";") if field.count("=")])
	return attrd

def get_gff_exons(gfffile, gffid, residue_offset=0):
	'''read GFF file, extract exon/CDS information for a specified gene, return a list of lists of residues as strings'''
	#TODO child_to_parent = {} # 

	# read through GFF, pull out all CDS features matching the ID, store length in bases
	linecounter = 0
	exoncounter = 0
	exon_lengths = {}
	print >> sys.stderr, "# Parsing GFF file {} for {}".format( gfffile, gffid ), time.asctime()
	for line in open(gfffile,'r'):
		line = line.strip()
		if line and line[0]!="#":
			linecounter += 1
			lsplits = line.split("\t")
			attrstring = lsplits[8]
			attr_dict = attributes_to_dict(attrstring)
			if attr_dict.get("ID",False)==gffid or attr_dict.get("Parent",False)==gffid:
				feature = lsplits[2]
				strand = lsplits[6]
				if feature == "CDS":
					exoncounter += 1
					exonlen = int(lsplits[4]) - int(lsplits[3]) + 1
					if strand=="+":
						exon_lengths[exoncounter] = exonlen
					elif strand=="-": # make numbers negative, so will be sorted backwards
						exon_lengths[-1*exoncounter] = exonlen
					else: # for anything else, as everything should have a strand
						print >> sys.stderr, "# WARNING UNKNOWN STRAND: {} FOR\n{}".format(strand, line)
	print >> sys.stderr, "# Finished parsing {}, read {} lines".format(gfffile, linecounter), time.asctime()
	if exoncounter==0:
		print >> sys.stderr, "# WARNING COULD NOT FIND ANY EXONS FOR {}".format(gffid)

	# sort CDS features in order, depending on strand, and count amino acids per exon
	basecounter = 0
	residuecounter = residue_offset
	cds_list = [] # list, of lists of residue numbers, as strings
	out_of_phase = [] # all residues that are out of phase, colored grey
	for exonnumber in sorted(exon_lengths.keys()):
		cdslength = exon_lengths[exonnumber]
		basecounter += cdslength
		firstresidue = residuecounter + 1
		residuecounter = (basecounter // 3) # add all integer residues
		# +1 for range, as last number is not included
		residue_strings = [str(i) for i in range(firstresidue,residuecounter+1)]
		cds_list.append(residue_strings)
		# if this exon is phase 1 or 2
		# and does not get back in phase from previous exons
		if (basecounter % 3):
			residuecounter += 1 # add 1 for the out of phase
			out_of_phase.append( str(residuecounter) )

	print >> sys.stderr, "# Found {} exons, {} bases for {} residues".format( exoncounter, basecounter, residuecounter )
	print >> sys.stderr, "# {} codons are out of phase".format( len(out_of_phase) )
	return cds_list, out_of_phase

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-c","--chain", default="A", help="chain used in PDB format file [A]")
	parser.add_argument("-g","--gff", help="GFF file that contains genes and CDS entries for -s", required=True)
	parser.add_argument("-i","--id", help="ID of sequence in the GFF file, if ID or Parent field does not match -s")
	args = parser.parse_args(argv)

	residue_list, outofphase_list = get_gff_exons(args.gff, args.id)
	make_output_script(residue_list, outofphase_list, args.chain, wayout)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
