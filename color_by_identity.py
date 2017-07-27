#!/usr/bin/env python
#
# color_by_identity.py v1 2017-07-25

'''color_by_identity.py v1 2017-07-26
adapted from consurf_new.py

in PyMOL, use by:

run color_by_identity.py
'''

def color_identity(selection="all", num_colors=9, MINIMUM=1.0, MAXIMUM=9.0):
	'''function to color atoms/residues by the temperatureFactor field'''

	# each color is RGB triplet, values ranging from 0 to 1, where black is 0 and white is 1
	# colors should be:
	# gray, light blue, blue, teal, green, yellow, orange, red, magenta
	# to correspond to identity of <= :
	# 0   , 50        , 60  , 70  , 80   , 90    , 95    , 98 , 100
	colors = [ [0.5, 0.5, 0.5] , [0.33,0.33,0.72] , [0,0.30,0.70], 
			  [0.06,0.68,0.40] , [0.21,0.66,0]    , [0.85,0.66,0], 
			  [0.94,0.32,0]    , [0.88,0,0]       , [1,0,0.55] ]

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
	cmd.select("insufficient", selection + " & b = 0")
	cmd.color("insufficient_color", "insufficient")

# This is required to make command available in PyMOL 
cmd.extend("color_identity", color_identity)

color_identity()
