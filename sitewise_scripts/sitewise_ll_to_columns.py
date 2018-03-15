#!/usr/bin/env python
#
# sitewise_ll_to_columns.py

'''
sitewise_ll_to_columns.py RAxML_perSiteLLs.Borowiec_best108_site_lk
'''

import sys
from collections import defaultdict

if len(sys.argv) < 2:
	print >> sys.stderr, __doc__
else:
	linecounter = 0
	sitelldict = defaultdict(list)
	headerline = ["site"]
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if linecounter < 1: # skip first line
			linecounter += 1
		else:
			header, sites = line.split("\t",1)
			headerline.append(header)
			for i,site in enumerate( sites.split(' ')):
				sitelldict[i].append(site)
	print >> sys.stdout, "\t".join(headerline)
	for k in sorted(sitelldict.keys()):
		print >> sys.stdout, "{}\t{}".format(k+1, "\t".join(sitelldict[k]))
