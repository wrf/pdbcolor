#!/usr/bin/env python
#
# color_by_likelihood.py v1 2018-02-01

'''color_by_likelihood.py v1 2018-02-01

in PyMOL, use by:

run color_by_likelihood.py
'''

def color_likelihood(selection="all", num_colors=18, MINIMUM=-1.0, MAXIMUM=16.0):
	'''function to color atoms/residues by the temperatureFactor field'''

	# each color is RGB triplet, values ranging from 0 to 1, where black is 0 and white is 1
	# colors should be:
	# for values [ -1 gray
	#            0 blue   1   2
	#            3 teal   4   5 sea-green
	#            6        7 green-brown
	#            8 rust   9 pink   a
	#            b rose   c        d
	#            e red    f
	#            16 green ]

	colors = [ [0.39, 0.39, 0.39] ,
             [0.00,0.22,0.72] , [0.05,0.38,0.71] , [0.07,0.47,0.70] ,
             [0.11,0.63,0.69] , [0.16,0.82,0.67] , [0.26,0.77,0.64] ,
             [0.39,0.72,0.60] , [0.42,0.49,0.42] ,
             [0.49,0.36,0.40] , [0.72,0.33,0.56] , [0.76,0.26,0.56] ,
             [0.78,0.24,0.48] , [0.77,0.22,0.44] , [0.77,0.20,0.33] ,
             [0.77,0.18,0.22] , [0.77,0.16,0.14] ,
             [0.21,0.66,0.00] ]

	bin_size = ((MAXIMUM - MINIMUM) + 1) / num_colors

	# iterate through temperatureFactor value in PyMOL
	# atoms with each value are selected and colored with the following steps
	for i in range(num_colors):

		lower = MINIMUM + i * bin_size
		upper = lower + bin_size - 0.01

		# Print out B-factor limits and the color for this group
		print lower, " - ", upper, " = ", colors[i]

		# Define a unique name for the atoms which fall into this group
		groupname = selection + "_group_" + str(i+1)

		# Compose a selection command which will select all atoms with some beta value
		sel_string = selection + " & ! b < " + str(lower)

		if(i < num_colors - 1):
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

	# color atoms with B-factor of 0 using the new color
	cmd.select("insufficient", selection + " & b > 11")
	cmd.color("insufficient_color", "insufficient")

# This is required to make command available in PyMOL 
cmd.extend("color_likelihood", color_likelihood)

color_likelihood()
