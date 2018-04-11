#!/usr/bin/env python
#
# color_by_identity.py v1 2017-07-25

'''color_by_identity.py v1.1 2018-04-02

    assumes that b-factor values have been recoded by pdb_site_identity.py

    in PyMOL, use by:

run color_by_identity.py

    to use a specific color palette, such as red:

color_by_identity(colorscheme='red')

    if multiple alignments are used for different proteins in a heteromultimer
    use the option bychain=True

color_by_identity(bychain=True)
'''

def color_by_identity(selection="all", colorscheme="rainbow", bychain=False):
	'''function to color atoms/residues by the temperatureFactor field'''

	# each color is RGB triplet, values ranging from 0 to 1, where black is 0 and white is 1

	colors = {} # key is colorscheme name, value is list of colors

	colors["red"] = [ [0.63,0.63,0.63] , [0.73,0.55,0.55] , [0.75,0.47,0.47], 
				   [0.77,0.38,0.38] , [0.79,0.29,0.29] , [0.82,0.21,0.21], 
				   [0.84,0.13,0.13]    , [0.88,0,0]     , [1,0,0.55] ]
	colors["green"] = [ [0.63,0.63,0.63] , [0.50,0.68,0.56] , [0.42,0.71,0.53], 
				   [0.35,0.74,0.49] , [0.26,0.77,0.44] , [0.19,0.80,0.41], 
				   [0.12,0.83,0.37] , [0.01,0.87,0.31] , [1,0,0.55] ]
	colors["blue"] = [ [0.63,0.63,0.63] , [0.50,0.58,0.68] , [0.42,0.55,0.71], 
				   [0.35,0.52,0.73] , [0.28,0.49,0.76] , [0.20,0.46,0.80], 
				   [0.12,0.43,0.83] , [0.00,0.38,0.87] , [1,0,0.55] ]
	colors["yellow"] = [ [0.63,0.63,0.63] , [0.66,0.67,0.51] , [0.68,0.71,0.43], 
				   [0.70,0.73,0.36] , [0.72,0.76,0.28] , [0.74,0.80,0.20], 
				   [0.76,0.83,0.12] , [0.79,0.87,0.00] , [1,0,0.55] ]
	# rainbow colors should be:
	# gray, light blue, blue, teal, green, yellow, orange, red, magenta
	# to correspond to identity of <= :
	# 0   , 50        , 60  , 70  , 80   , 90    , 95    , 98 , 100
	colors["rainbow"] = [ [0.5, 0.5, 0.5] , [0.38,0.38,0.72] , [0,0.30,0.70], 
			  [0.06,0.68,0.40] , [0.21,0.66,0]    , [0.85,0.66,0], 
			  [0.94,0.32,0]    , [0.88,0,0]       , [1,0,0.55] ]

	groupnames = [ "0", "50pct", "60pct", "70pct", "80pct", "90pct", "95pct", "98pct", "100pct" ]

	num_colors = 9
	# need 101 as extra unused value to prevent IndexError
	binvalues = [00.0 ,50.0 ,60.0 ,70.0 ,80.0 ,90.0 ,95.0 ,98.0 ,100, 101]

	# separately color individual chains when multiple alignments are used
	if bychain:
		schemelist = ["red","green","blue","yellow"]
		chain_to_obj = {}
		chain_to_colorscheme = {}
		for obj in cmd.get_object_list("all"):
			chainlist = cmd.get_chains(obj)
			for i,chain in enumerate(chainlist):
				chain_to_obj[chain] = "chain "+chain
				chain_to_colorscheme[chain] = schemelist[(i+4)%4]
				print "chain " + chain + " as colorscheme " + schemelist[(i+4)%4]
	else: # color all chains, or the only chain
		chain_to_obj = {selection:selection}
		if colorscheme not in colors: # default is rainbow
			print colorscheme + " is not an option, using 'rainbow' scheme"
			chain_to_colorscheme = {selection: "rainbow" }
		else: # assign scheme given by user
			chain_to_colorscheme = {selection: colorscheme }

	# iterate through temperatureFactor value in PyMOL
	# atoms with each value are selected and colored with the following steps
	for chain in chain_to_colorscheme.iterkeys():
		for i in range(num_colors):

			lower = binvalues[i]
			upper = binvalues[i+1]-0.0001

			# Define a unique name for the atoms which fall into this group
			groupname = groupnames[i] + "_grp_" + str(i+1) + chain

			# Compose a selection command which will select all atoms with some beta value
			sel_string = chain_to_obj[chain] + " & ! b < " + str(lower)

			# Print out B-factor limits and the color for this group
			if(i < num_colors - 1):
				print lower, " - ", upper, " = ", colors[chain_to_colorscheme[chain]][i]
				sel_string += " & b < " + str(upper)
			else:
				print " >= ", lower, " = ", colors[chain_to_colorscheme[chain]][i]
				#sel_string += " & ! b > " + str(upper)

			# select the atoms
			cmd.select(groupname, sel_string) 

			# create a new color where name is color_number and values come from the colors list
			color_name = "color_" + str(i+1) + chain
			cmd.set_color(color_name, colors[chain_to_colorscheme[chain]][i])

			# this step actually recolors the residues
			# all selected atoms are reassigned the color based on the name from the step above
			cmd.color(color_name, groupname)

	# this is an extra color for errors, etc, leftover from consurf_new.py
	# Create new color for insufficient sequences
	insuf_color = [0.75, 0.75, 0.58823529]
	cmd.set_color("insufficient_color", insuf_color)

	# color atoms with B-factor of 0 using the new color
	cmd.select("insufficient", selection + " & b < 0")
	cmd.color("insufficient_color", "insufficient")

# This is required to make command available in PyMOL 
cmd.extend("color_by_identity", color_by_identity)

color_by_identity()

print "colorscheme can be changed to 'red' 'green' 'blue' or yellow'"
print "this is called as:  color_by_identity(colorscheme='red')"
print "to color multiple chains separately, use:  color_by_identity(bychain=True)"
