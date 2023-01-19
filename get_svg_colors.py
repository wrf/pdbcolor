#!/usr/bin/env python
#
# get_svg_colors.py 2019-09-25

'''get_svg_colors.py  last modified 2023-01-19

  this script was mostly intended to extract colors from a figure
  used to generate a color scheme in PyMOL

  it will extract all fill colors used in an .svg file
  and print the hex values converted to 0-1 scale for PyMOL coloring

get_svg_colors.py image.svg

  output is printed as:
62698d 0.38 0.41 0.55
  for hex color, then 0.0-1.0 values of R G and B

  used in PyMOL like:
set_color green98, [0.01,0.87,0.31]

'''

import re
import sys

if len(sys.argv)<2:
	sys.stderr.write( __doc__ )
else:
	for line in open(sys.argv[1]):
		fillcolor = re.search("fill:#(\w{6});", line)
		if fillcolor:
			hexcolor = fillcolor.group(1)
			sys.stdout.write( "{}  {:.2f},{:.2f},{:.2f}\n".format(hexcolor, *[ int(color,16)/255.0 for color in [ hexcolor[0:2] , hexcolor[2:4] , hexcolor[4:6] ] ]) )
