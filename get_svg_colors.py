#!/usr/bin/env python
#
# get_svg_colors.py

'''
  extract all fill colors used in an svg

get_svg_colors.py image.svg

  output is printed as:
62698d 0.38 0.41 0.55
  for hex color, then 0.0-1.0 values of R G and B
'''

import re
import sys

if len(sys.argv)<2:
	print >> sys.stderr, __doc__
else:
	for line in open(sys.argv[1]):
		fillcolor = re.search("fill:#(\w{6});", line)
		if fillcolor:
			hexcolor = fillcolor.group(1)
			print >> sys.stdout, hexcolor, "{:.2f},{:.2f},{:.2f}".format(*[ int(color,16)/255.0 for color in [ hexcolor[0:2] , hexcolor[2:4] , hexcolor[4:6] ] ])
