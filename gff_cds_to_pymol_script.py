#!/usr/bin/env python
#
# gff_cds_to_pymol_script.py v1 2018-06-04

'''gff_cds_to_pymol_script.py  last modified 2022-02-01
    script to generate a PyMOL script to color residues corresponding to exons
    20 colors are currently available, and repeat after 20 exons

gff_cds_to_pymol_script.py -g ppyr_00001.gff -i PPYR_00001-PA > ppyr_00001_color_5dv9.pml

    the script detects CDS features in a GFF where the ID matches -i
    this assumes the first codon in the CDS feature maps to residue 1
    even if those residues from that exon are not in the structure

    exons are automatically numbered, corresponding to the first CDS feature
    this might not match the first exon in the gene, due to alternative splicing
    or untranslated exons

    default chain is A, and can be changed with -c
    for homomultimers, specify the chains as -c AB

gff_cds_to_pymol_script.py -g human_IDH2.gff -i cds82740 -c AB > IDH2_color_5i96.pml

run in PyMOL as (meaning do not use the run command):
@ppyr_00001_color_5dv9.pml
'''

import sys
import os
import time
import argparse

def make_output_script(cds_to_residues, out_of_phase_list, chainset, wayout):
	'''from the GFF information, print a script for PyMOL'''

	# colors in order of:
	# blue, green, red, orange, purple
	# brown, teal, lilac, seagreen, yellow
	exon_colors = { 1 :[0.65,0.81,0.89], 2 :[0.12,0.47,0.71], 
                    3 :[0.70,0.87,0.54], 4 :[0.20,0.63,0.17], 
                    5 :[0.98,0.60,0.60], 6 :[0.89,0.10,0.11], 
                    7 :[0.99,0.75,0.44], 8 :[1.00,0.50,0.00], 
                    9 :[0.79,0.70,0.84], 10 :[0.42,0.24,0.60], 
                    11 :[0.90,0.65,0.50], 12 :[0.67,0.35,0.14], 
                    13 :[0.19,0.91,0.94], 14 :[0.00,0.72,0.75], 
                    15 :[0.89,0.65,0.82], 16 :[0.71,0.12,0.48], 
                    17 :[0.65,0.89,0.76], 18 :[0.12,0.71,0.45], 
                    19 :[0.89,0.88,0.65], 20 :[0.71,0.67,0.12] }
	# { 1:"#a6cee3", 2:"#1f78b4", 3:"#b2df8a", 4:"#33a02c", 5:"#fb9a99", 
	# 6:"#e31a1c", 7:"#fdbf6f", 8:"#ff7f00", 9:"#cab2d6", 10:"#6a3d9a", 
	# 11:"#e5a580", 12:"#aa5823", 13:"#30e7ef", 14:"#00b7c0", 15:"#e3a6d1", 
	# 16:"#b41f7a", 17:"#a6e3c1", 18:"#1fb472", 19:"#e3e0a6", 20:"#b4ac1f" }

	outofphase_color = [0.3,0.3,0.3] # dark grey
	numcolors = len(exon_colors)
	chainstring = " or ".join( ["chain {}".format(chain) for chain in chainset] )

	# begin printing commands for PyMOL script
	wayout.write("hide everything\n")
	wayout.write("bg white\n")
	wayout.write("show cartoon, ({})\n".format(chainstring) )
	for exonnum, residues in enumerate(cds_to_residues):
		colorlist = map(str, exon_colors[(exonnum % numcolors)+1])
		wayout.write("set_color excolor_{}, [{}]\n".format( exonnum+1, ",".join(colorlist) ) )
		wayout.write("select exon_{}, (resi {}) & ({})\n".format( exonnum+1, "-".join(residues), chainstring ) )
		wayout.write("color excolor_{}, exon_{}\n".format( exonnum+1, exonnum+1 ) )
	outofphase_col_str = map(str, outofphase_color)
	wayout.write("set_color outphase, [{}]\n".format( ",".join(outofphase_col_str) ) )
	wayout.write("select out_of_phase, (resi {}) & ({})\n".format( ",".join(out_of_phase_list), chainstring ) )
	wayout.write("color outphase, out_of_phase\n")
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
	cdscounter = 0
	exon_lengths = {}
	gff_ids_from_file = {} # list of the IDs found in the file, for error checking, key is ID, value is True
	sys.stderr.write("# Parsing GFF file {} for {}  {}\n".format( gfffile, gffid , time.asctime() ) )
	for line in open(gfffile,'r'):
		line = line.strip()
		if line and line[0]!="#":
			linecounter += 1
			lsplits = line.split("\t")
			attrstring = lsplits[8]
			attr_dict = attributes_to_dict(attrstring)
			cds_id = attr_dict.get("ID",False)
			if cds_id:
				gff_ids_from_file[cds_id] = True
			if cds_id==gffid or attr_dict.get("Parent",False)==gffid:
				feature = lsplits[2]
				strand = lsplits[6]

				if feature == "CDS":
					cdscounter += 1
					exonlen = int(lsplits[4]) - int(lsplits[3]) + 1
					if strand=="+":
						exon_lengths[cdscounter] = exonlen
					elif strand=="-": # make numbers negative, so will be sorted backwards
						exon_lengths[-1*cdscounter] = exonlen
					else: # for anything else, as everything should have a strand
						sys.stderr.write("# WARNING UNKNOWN STRAND: {} FOR\n{}\n".format(strand, line) )
	sys.stderr.write("# Finished parsing {}, read {} lines  {}\n".format(gfffile, linecounter, time.asctime() ) )
	if cdscounter==0: # could not find any CDS features, or ID is wrong
		sys.stderr.write("# WARNING COULD NOT FIND ANY CDS FEATURES FOR {}\n".format(gffid) )
		if gff_ids_from_file: # list IDs that were found
			sys.stderr.write("# FILE CONTAINS ID={}\n".format( " ".join(gff_ids_from_file.keys()) ) )

	# sort CDS features in order, depending on strand, and count amino acids per exon
	basecounter = 0
	residuecounter = residue_offset
	cds_list = [] # list, of lists of residue numbers, as strings
	out_of_phase = [] # all residues that are out of phase, colored grey
	# append a residue start-end pair to cds_list
	for exonnumber in sorted(exon_lengths.keys()):
		cdslength = exon_lengths[exonnumber]
		basecounter += cdslength
		firstresidue = residuecounter + 1
		residuecounter = (basecounter // 3) # add all integer residues
		# make a list of the first and last residues, to become strings
		residue_strings = map(str, [firstresidue, residuecounter] )
		# append the 2-item-list to cds_list
		cds_list.append(residue_strings)
		# if this exon is phase 1 or 2
		# and does not get back in phase from previous exons
		if (basecounter % 3):
			residuecounter += 1 # add 1 for the out of phase
			out_of_phase.append( str(residuecounter) )

	if cdscounter:
		sys.stderr.write("# Found {} CDS exons, {} bases for {} residues\n".format( cdscounter, basecounter, residuecounter ) )
		sys.stderr.write("# {} codons are out of phase\n".format( len(out_of_phase) ) )
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
