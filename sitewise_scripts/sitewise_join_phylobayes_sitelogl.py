#!/usr/bin/env python
#
# sitewise_join_phylobayes_sitelogl.py created by WRF 2018-03-12

'''
sitewise_join_phylobayes_sitelogl.py tree1-chain1.sitelogl tree2-chain1.sitelogl tree3-chain1.sitelogl > combined_sitelogl.tab

    generate sitelogl with:
mpirun -np 6 ~/phylobayes/pbmpi/data/readpb_mpi -sitelogl -x 50 5 tree1-chain1-1
'''

import sys
from collections import defaultdict

if len(sys.argv) < 2:
	print >> sys.stderr, __doc__
else:
	linecounter = 0
	loglsbysite = defaultdict(list)
	siteloglfiles = sys.argv[1:]
	header = "site\t{}".format( "\t".join( ["T{}".format(i+1) for i in range(len(siteloglfiles)) ] ) )
	for siteloglfile in siteloglfiles:
		print >> sys.stderr, "# reading sitelogl from {}".format(siteloglfile)
		for line in open(siteloglfile, 'r'):
			line = line.strip()
			if line:
				lsplits = line.split("\t")
				position = int(lsplits[0])
				sitelogl = lsplits[1]
				CPOindex = lsplits[2] # conditional predictive ordinate from Lewis 2013 Syst Biol
				loglsbysite[position].append(sitelogl)
	print >> sys.stdout, header
	for k in sorted(loglsbysite.keys()):
		print >> sys.stdout, "{}\t{}".format( k, "\t".join(loglsbysite[k]) )
