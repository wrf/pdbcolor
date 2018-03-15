#!/usr/bin/env python

# sitewise_columns_to_fasta.py

'''sitewise_columns_to_fasta.py last modified 2018-03-15
    convert tabular sitewise log-likelihood to a FASTA string

sitewise_columns_to_fasta.py RAxML_perSiteLLs.matrix1.tab > RAxML_perSiteLLs.matrix1.fasta

    by default, this assumes two or three trees (called T1, T2, and T3)
    difference in lnL is ln(max)-ln(median),
    for whichever trees are max and median, respectively
    values are then coded to a single character for the ranges:
    0-0.5 as 0, 0.5-2.0 as 1, and >2.0 as 2 for each tree

    values 0-2 favor T1, 3-5 favor T2, and 6-8 favor T3

    change the strong site cutoff (0.5) with an additional term:

sitewise_columns_to_fasta.py RAxML_perSiteLLs.matrix1.tab 1.0

    generate the tabular sitewise results using:
sitewise_ll_to_columns.py RAxML_perSiteLLs.matrix1 > RAxML_perSiteLLs.matrix1.tab
'''

import sys
from collections import defaultdict
from numpy import median

if len(sys.argv) < 2:
	print >> sys.stderr, __doc__
else:
	treecounter = defaultdict(int)
	strongcounter = defaultdict(int)
	constsites = 0
	totalsites = 0
	fastastring = ""
	strongsitecutoff = 0.5
	if len(sys.argv) > 2: # assume alternate strong cutoff is used
		print >> sys.stderr, "# Using {} as strong site cutoff".format(sys.argv[2])
		strongsitecutoff = float(sys.argv[2])
	print >> sys.stderr, "# Converting {} to a fasta string".format(sys.argv[1])
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			if lsplits[0] == "site":
				continue
			totalsites += 1
			if lsplits[1] == "const":
				fastastring += "x"
				constsites += 1
			else:
				lnlfloats = [float(n) for n in lsplits[1:]]
				# was calculated before as mean of T1-T2, T2-T3, T1-T3
				# but should instead be max - median
				# to avoid cases where min is much worse than max and med
				dlnl = max(lnlfloats) - median(lnlfloats) # should always be greater than 0
				toptree = lnlfloats.index(max(lnlfloats))
				treecounter[toptree+1] += 1
				if dlnl < strongsitecutoff: # should be 0.5 by default
					adjdlnl = 0
				else: # meaning above noise threshold
					strongcounter[toptree+1] += 1
					if dlnl < 2.0: # strong site, not very strong
						adjdlnl = 1
					else: # very strong sites
						adjdlnl = 2
				adjdlnl = adjdlnl + (3 * toptree)
				fastastring += str(adjdlnl)
	print >> sys.stdout, ">diffLnL\n{}".format(fastastring)
	print >> sys.stderr, "# {} total sites".format(totalsites)
	print >> sys.stderr, "# {} constant sites".format(constsites)
	for k in sorted(treecounter.keys()):
		print >> sys.stderr, "# {} total sites, {} strong sites, for T{}".format(treecounter[k], strongcounter[k], k)
