#!/usr/bin/env python
#
# pdb_color_generic.py v1 2019-01-21

'''pdb_color_generic.py  last modified 2019-02-14
    generate a script to color a PDB file based on generic tabular data
    REQUIRES numpy for arange

pdb_color_generic.py -c 4 -d , -p 3rze.pdb -i 3rze.map.rates_features.csv -l green > 3rze.color_by_dnds.pml

    options are:
    -c column, column index of desired data, starting from 0
    -d delimiter, default is tab, e.g. change to "," for csv
    -g group name, changes the group name in the PyMOL selection
    -l color schemes, schemes include:
      sequential - from gray to intense color (e.g. -l green )
        red, yellow, blue, green
      diverging - with either light or dark in the middle ( -l div2b )
      so values close to 0 will blend in with the background color
        div1w - dark brown to white to green, for white bg display
        div1b - light brown to black to light green, for black bg display
        div2w - brick red to white to dark blue, for white bg
        div2b - pink to black to sky blue, for black bg

    data -i can be any text file, so long as columns can be split with -d
    generally, it is best to format in something like this:
residue score   chain
1       88.5    A

    chain information is optional, and the residue column can be specified
    with the option --site-column, starting from 0

    if multiple proteins are in the PDB file, then it is better to have
    the chain specified in the data file. Otherwise the script can be run
    multiple times, each time specifying a different chain by --default-chain
'''

import sys
import argparse
import itertools
import numpy # for arange,around,floor,log10,median
from collections import defaultdict

def read_generic_data(datafile, delimiter, scorecolumn, sitecolumn=0, chaincolumn=None):
	'''read generic data file, and return a dict key is site number and value is score'''
	rawscore_dict = defaultdict(dict) # key is chain, value is dict of site and raw score
	linecounter = 0
	foundscores = 0
	print >> sys.stderr, "# Reading scores from column {} in {}, separated by {}".format( scorecolumn, datafile, delimiter )
	for line in open(datafile, 'r'):
		line = line.strip()
		if line and line[0]!="#":
			linecounter += 1
			lsplits = line.split(delimiter)
			sitecol = lsplits[sitecolumn] # raw, might be string
			try: # check if can be turned into integer
				site = int(float(sitecol))
			except ValueError: # no value, row is empty in this column or is header
				continue
			score = float(lsplits[scorecolumn])
			if type(chaincolumn) is int:
				chain = lsplits[int(chaincolumn)]
			else:
				chain = "A"
			rawscore_dict[chain][site] = float(score)
			foundscores += 1
	print >> sys.stderr, "# Counted {} lines with {} scores".format( linecounter, foundscores )
	return rawscore_dict

def get_chains_only(defaultchain, seqidlist, pdbfile):
	'''read PDB file and return a dict where key is chain and value is sequence ID'''
	keepchains = {} # dict where key is chain and value is seqid
	print >> sys.stderr, "# Reading chain from PDB {}".format(pdbfile)
	for line in open(pdbfile,'r'):
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
	if defaultchain: # meaning nothing was found, use default and single sequence
		if seqidlist:
			keepchains[defaultchain] = seqidlist[0]
		else:
			keepchains[defaultchain] = "UNKNOWN"
		print >> sys.stderr, "### using default chain {}".format( defaultchain )
	return keepchains

def make_output_script(wayout, scoredict, keepchains, groupname, exclude_first, basecolor="red", reverse_colors=False, ZEROOVERRIDE=0.0):
	'''from the identity calculations, print a script for PyMOL'''
	###
	### DECLARE COLOR SCHEMES ###
	###
	colors = {} # key is colorscheme name, value is list of colors
	# sequential color schemes
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

	# diverging scheme 1, where neutral is white
	# brown to white to dark green
	#543005 #8c510a #bf812d #dfc27d #f6e8c3
	#f5f5f5
	#c7eae5 #80cdc1 #35978f #01665e #003c30
	colors["div1w"] = [ [0.55,0.32,0.04] , [0.75,0.51,0.18] , [0.87,0.76,0.49] ,
				   [0.96,0.91,0.76] , [0.96,0.96,0.96] , [0.78,0.92,0.90] , 
				   [0.50,0.80,0.76] , [0.21,0.59,0.56] , [0.00,0.40,0.37] ]

	# diverging scheme 1, where neutral is black
	colors["div1b"] = [ [0.87,0.76,0.49] , [0.75,0.51,0.18] , [0.55,0.32,0.04] , 
				   [0.33,0.19,0.02] , [0.24,0.24,0.24] , [0.00,0.40,0.37] ,
				   [0.21,0.59,0.56] , [0.50,0.80,0.76] , [0.78,0.92,0.90] ]

	# diverging scheme 2, neutral is white
	# brick red to white to dark blue
	#67001f #b2182b #d6604d #f4a582 #fddbc7
	#f7f7f7
	#d1e5f0 #92c5de #4393c3 #2166ac #053061
	colors["div2w"] = [ [0.70,0.09,0.17] , [0.84,0.38,0.30] , [0.96,0.65,0.51] , 
				   [0.99,0.86,0.78] , [0.97,0.97,0.97] , [0.82,0.90,0.94] , 
				   [0.57,0.77,0.87] , [0.26,0.58,0.76] , [0.13,0.40,0.67] ]

	# diverging scheme 2, neutral is black for white bg
	# pink to black to light blue
	colors["div2b"] = [ [0.96,0.65,0.51] , [0.84,0.38,0.30] , [0.70,0.09,0.17] , 
				   [0.40,0.00,0.12] , [0.24,0.24,0.24] , [0.02,0.19,0.38] , 
				   [0.13,0.40,0.67] , [0.26,0.58,0.76] , [0.57,0.77,0.87] ]

	insuf_color = [0.75, 0.75, 0.58]

	###
	### DETERMINE OPTIMAL BIN VALUES DIRECTLY FROM DATA ###
	###
	all_scores = list( itertools.chain( rd.values() for rd in scoredict.values() ) )[0]
	print >> sys.stderr, "# Generating list of bins from data"
	lowest_val = min(all_scores)
	highest_val = max(all_scores)
	val_range = highest_val - lowest_val
	magnitude = int(numpy.floor(numpy.log10(val_range)))

	round_range = numpy.around([lowest_val,highest_val], decimals=abs( min([magnitude,0])) )
	rounded_diff = round_range[1]-round_range[0]

	median_val = numpy.median(all_scores)

	#print >> sys.stderr, magnitude, round_range, rounded_diff
	print >> sys.stderr, "# data range from {:.2f} to {:.2f}, diff of {:.2f}, median of {:.2f}".format(lowest_val, highest_val, val_range, median_val)

	# correction if rounded lower bound is greater than the lowest value
	round_correction = 10**magnitude
	if round_range[0] > lowest_val:
		print >> sys.stderr, "### correcting lower bound {:.2f} to {:.2f}".format( round_range[0], round_range[0]-round_correction)
		round_range[0] = round_range[0] - round_correction

	if 1 < val_range < 8:
		magnitude = magnitude-1
		print >> sys.stderr, "### value range is lower than 8, adjusting decimals in selection names from {} to {}".format(magnitude, magnitude-1)

	# all values are positive
	if lowest_val >= 0:
		print >> sys.stderr, "### lowest value is 0, best use sequential colors"
		last_bin = float(round_range[1] + round_correction)
		rounded_step = rounded_diff/8
		binvalues = numpy.arange(round_range[0], round_range[1]+rounded_step, rounded_step).tolist()
		binvalues.append( last_bin )
		if median_val < binvalues[1] and exclude_first is False:
			print >> sys.stderr, "### median value is low: {:.2f} , use --exclude-first-group if there are too many residues in group 1".format(median_val)

	# values span 0, set up so zero is middle
	# ZEROOVERRIDE is 0 by default
	elif lowest_val < 0 and highest_val > 0:
		print >> sys.stderr, "### values bridge 0, best use diverging color schemes"
		# middle value should be 0, thus make two ranges
		low_step = abs(round_range[0]/4)
		# last value in list should be ZEROOVERRIDE, 0
		low_set = numpy.arange(round_range[0], ZEROOVERRIDE+round_correction/10, low_step).tolist()
		high_step = round_range[1]/4
		# first value should be slightly above zero, to catch near zero values
		high_set = numpy.arange(ZEROOVERRIDE+round_correction/10, round_range[1]+high_step, high_step).tolist()
		#print >> sys.stderr, low_set, low_step, high_set, high_step
		binvalues = low_set + high_set
	elif highest_val <= 0: # all values are negative
		print >> sys.stderr, "### all values are negative, best use sequential colors and --reverse-colors"

	# make correction factor for naming to avoid decimals
	binname_correction = 10/10**magnitude

	###
	### PRINT COMMANDS FOR SCRIPT ###
	###
	print >> sys.stderr, "# Generating PyMOL script for color scheme {} with bins of:\n{}".format( basecolor, binvalues )
	print >> wayout, "hide everything"
	#print >> wayout, "bg white"
	print >> wayout, "show cartoon"
	print >> wayout, "set_color colordefault, [{}]".format( ",".join(map(str,insuf_color)) )
	print >> wayout, "color colordefault, all"
	# make commands for target color
	targetcolors = colors[basecolor]
	if reverse_colors:
		targetcolors.reverse()
	for i,rgb in enumerate(targetcolors):
		colorname = "{}{:02d}".format( basecolor, int(binvalues[i]*binname_correction) )
		print >> wayout, "set_color {}, [{}]".format( colorname, ",".join(map(str,rgb)) )

	# make commands for each chain
	for chain in keepchains.iterkeys(): # keys are chain letters, values are seq IDs
		scoregroups = defaultdict(list) # key is percent group, value is list of residues
		# for each residue, assign to a bin
		for residue in scoredict[chain].iterkeys():
			residuescore = scoredict[chain].get(residue,0.00)
			for i,value in enumerate(binvalues[:-1]):
				upper = binvalues[i+1]
				if residuescore < upper:
					scoregroups[value].append(residue)
					break
			# should not need an else if last bin is large enough
		# assign whole chain to lowest color, then build up
		print >> wayout, "color {}{:02d}, chain {}".format( basecolor, int(value*binname_correction), chain )
		# for each bin, make a command to color all residues of that bin
		for i,value in enumerate(binvalues[:-1]):
			if i==0 and exclude_first: # long lists apparently crash the program, so skip
				continue
			binname = "{:02d}_{}_{}_{}".format( int(value*binname_correction), groupname, i+1, chain )
			resilist = map(str,scoregroups[value])
			if resilist: # do not print empty groups
				binresidues = ",".join(resilist)
				print >> wayout, "select {}, (chain {} & resi {})".format( binname, chain, binresidues )
				print >> wayout, "color {}{:02d}, {}".format( basecolor, int(value*binname_correction), binname )
	# no return

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-c","--data-column", default=1, type=int, help="index of data column, starting from 0 [1]")
	parser.add_argument("-d","--delimiter", default="\t", help="delimiter for file, default is tab")
	parser.add_argument("-g","--group-name", default="grp", help="name for groups, default is grp, appears as 100_grp_9_A")
	parser.add_argument("-i","--input-file", help="tabular or csv file of site data, with sites in the first column", required=True)
	parser.add_argument("-l","--base-color", default="red", help="color gradient, default is red, options are: [red,yellow,blue,green,div1w,div1b,div2w,div2b]")
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	parser.add_argument("-s","--sequence", nargs="*", help="sequence ID for PDB, give multiple names if data is available in the input file")
	parser.add_argument("--default-chain", default="A", help="default letter of chain [A], if DBREF for the sequence cannot be found in PDB")
	parser.add_argument("--exclude-first-group", action="store_true", help="exclude lowest score group, for cases where there are a large number of score-0 residues")
	parser.add_argument("--exclude-middle-group", action="store_true", help="exclude middle score group for diverging-type datasets")
	parser.add_argument("--reverse-colors", action="store_true", help="if used, reverse colors for negative-value datasets")
	parser.add_argument("--chain-column", help="column containing chain ID, default is None, will use chain A")
	parser.add_argument("--site-column", default=0, type=int, help="index of site column, starting from 0 [0]")
	parser.add_argument("--zero-override", default=0.0, type=float, help="middle index if data spans negative to positive, default is 0 as the middle color")
	args = parser.parse_args(argv)

	# read generic format data
	datadict = read_generic_data(args.input_file, args.delimiter, args.data_column, args.site_column, args.chain_column)

	# make PyMOL script with color commands
	refchains = get_chains_only(args.default_chain, args.sequence, args.pdb)
	make_output_script(wayout, datadict, refchains, args.group_name, args.exclude_first_group, args.base_color, args.reverse_colors, args.zero_override)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
