#!/usr/bin/env python
#
# color_by_likelihood.py v1 2018-02-01

'''color_by_likelihood.py v1.2 2022-02-04

in PyMOL, use by:

run color_by_likelihood.py
'''

def color_likelihood(selection="all", whitebg=False):
	'''function to color atoms/residues by the temperatureFactor field'''

	# each color is RGB triplet, values ranging from 0 to 1, where black is 0 and white is 1
	# colors should be:
	# for values [ -1 gray
	#            0 clay   1 pink   2 red
	#            3 mud    4 teal   5 green
	#            6 lead   7 cobalt 8 blue
	#            9 orange ]

	colors = [ [0.39, 0.39, 0.39] ,
             [0.49,0.36,0.40] , [0.76,0.26,0.56] , [0.77,0.20,0.33] ,
             [0.42,0.49,0.42] , [0.26,0.77,0.64] , [0.19,0.74,0.34] ,
             [0.39,0.37,0.55] , [0.35,0.26,0.76] , [0.24,0.45,0.87] ,
             [1.00,0.60,0.12] ]

	if whitebg: # replaces colors -1, 0, 3, 6, 9
		colors = [ [0.99,0.87,0.37] ,
             [0.83,0.74,0.76] , [0.76,0.26,0.56] , [0.77,0.20,0.33] ,
             [0.73,0.82,0.73] , [0.26,0.77,0.64] , [0.19,0.74,0.34] ,
             [0.74,0.73,0.86] , [0.35,0.26,0.76] , [0.24,0.45,0.87] ,
             [0.84,0.46,0.00] ]

	# colors for:
	# 7d7d00 dark yellow    b28b71 sand    640064 burgundy
	# 966400 orange    785583 lilac    4b6419 brown
	gapcolors = [ [0.49,0.49,0.00] , [0.70,0.55,0.44] , [0.39,0.00,0.39] ,
                 [0.59,0.39,0.00] , [0.47,0.33,0.51] , [0.29,0.39,0.10] ] 

	groupnames = [ "gaps", 
                  "t1-weak", "t1-strong", "t1-max",
                  "t2-weak", "t2-strong", "t2-max",
                  "t3-weak", "t3-strong", "t3-max",
                   "const" ]

	bin_size = 1
	NUMCOLORS = 11
	# iterate through temperatureFactor value in PyMOL
	# atoms with each value are selected and colored with the following steps
	for i in range(NUMCOLORS):

		lower = -1 + i
		upper = lower + 0.99

		# Print out B-factor limits and the color for this group
		print("{} = {}".format( lower, colors[i] ) )

		# Define a unique name for the atoms which fall into this group
		groupname = groupnames[i] + "_grp" + str(i+1)

		# Compose a selection command which will select all atoms with some beta value
		sel_string = selection + " & ! b < " + str(lower)

		if(i < NUMCOLORS - 1):
			sel_string += " & b < " + str(upper)
		else:
			sel_string += " & ! b > " + str(upper)

		# select the atoms
		cmd.select(groupname, sel_string)

		# create a new color where name is color_number and values come from the colors list
		color_name = "color_" + str(i+1)
		cmd.set_color(color_name, colors[i])

		# this step actually recolors the residues
		# all selected atoms are reassigned the color based on the name from the step above
		cmd.color(color_name, groupname)

	# this is an extra color for errors, etc, leftover from consurf_new.py
	# Create new color for insufficient sequences
	insuf_color = [0.75, 0.75, 0.58823529]
	cmd.set_color("insufficient_color", insuf_color)

	# color atoms with B-factor of 99 using the new color
	cmd.select("insufficient", selection + " & b > 9")
	cmd.color("insufficient_color", "insufficient")

	if whitebg:
		print("using color set for white backgrounds")
	else:
		print("to use the color scheme for white backgrounds, type: color_likelihood(whitebg=True)")

# This is required to make command available in PyMOL 
cmd.extend("color_likelihood", color_likelihood)

color_likelihood()
