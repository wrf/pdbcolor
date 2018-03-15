#!/usr/bin/env python
#
# alignment_pair_stats.py created 2017-10-17

'''alignment_pair_stats.py  last modified 2017-11-13

alignment_pair_stats.py -a align_pairs/

    alignment pairs can be generated using blast_to_align_pairs.py
'''

import sys
import os
import argparse
import time
import re
from glob import glob
from Bio import AlignIO

def check_alignment_pair(alignfile, alignformat):
	'''read alignment and extract information regarding coverage'''
	print >> sys.stderr, "# Reading alignment from {}".format( alignfile )
	alignment = AlignIO.read( alignfile, alignformat )
	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( num_taxa, al_length )
	trimmedlength = len( str(alignment[0].seq).replace("-","") )
	trimmedspan = re.search("\w[\w-]+\w", str(alignment[0].seq) ).span()
	spanlength = int(trimmedspan[1]) - int(trimmedspan[0])
	refprotlength = len( alignment[1].seq )
	statstring = "{}\t{}\t{}\t{:.3f}\t{}\t{}\t{:.3f}\t{}".format( alignment[0].id , alignment[1].id , trimmedlength , 1.0*trimmedlength/refprotlength, trimmedspan , spanlength, 1.0*spanlength/refprotlength, refprotlength )
	return statstring, trimmedlength, refprotlength

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignments', nargs="*", help="alignment files or directory")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-s','--sort', action="store_true", help="sort filenames before processing")
	args = parser.parse_args(argv)

	if os.path.isdir(args.alignments[0]):
		print >> sys.stderr, "# Reading alignments files from directory {}".format(args.alignments[0]), time.asctime()
		globstring = "{}*.aln".format(args.alignments[0])
		alignmentfiles = glob(globstring)
	elif os.path.isfile(args.alignments[0]):
		alignmentfiles = args.alignments
	else:
		raise OSError("ERROR: Unknown alignments, exiting")
	if args.sort:
		alignmentfiles = sorted(alignmentfiles, key=lambda x: int(re.search("(\d+)-(\d+)-",x).group(1)))

	headerline = "partition\tprotID\ttrimmedLength\ttrimmedPercent\tspan\tspanLength\tspanPercent\trefProtLength"
	print >> wayout, headerline
	refsitesum = 0
	trimsitesum = 0
	for alignfile in alignmentfiles:
		alignstats, trimsites, refsites = check_alignment_pair(alignfile, args.format)
		refsitesum += refsites
		trimsitesum += trimsites
		print >> wayout, alignstats
	print >> sys.stderr, "# Overall coverage is {} of {}, {:.3f} of reference sites".format( trimsitesum, refsitesum, 1.0*trimsitesum/refsitesum ), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
