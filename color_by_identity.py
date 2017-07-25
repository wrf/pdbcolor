#!/usr/bin/env python
#
# color_by_identity.py v1 2017-07-25

'''color_by_identity.py v1 2017-07-25
adapted from consurf_new.py

in PyMOL, use by:

run color_by_identity.py
'''

def color_identity(selection="all"):
	'''function to color atoms/residues by the temperatureFactor field'''
	MINIMUM = 1.0
	MAXIMUM = 9.0

	num_colors = 9
	# colors should be:
	# gray, light blue, blue, teal, green, yellow, orange, red, magenta
	# to correspond to identity of <= :
	# 0   , 50        , 60  , 70  , 80   , 90    , 95    , 98 , 100
	colors = [ 
	[0.5, 0.5, 0.5] ,
	[0.33,0.33,0.78], 
	[0,0.30,0.70], 
	[0.06,0.68,0.33], 
	[0.21,0.66,0], 
	[0.85,0.66,0], 
	[0.94,0.24,0], 
	[0.88,0,0],
	[1,0,0.55] ]
	
	bin_size = ((MAXIMUM - MINIMUM) + 1) / num_colors

	for i in range(num_colors):

		lower = MINIMUM + i * bin_size  
		upper = lower + bin_size - 0.01   
		color = colors[i]

		# Print out B-factor limits and the color for this group
		print lower, " - ", upper, " = ", color

		# Define a unique name for the atoms which fall into this group
		group = selection + "_group_" + str(i+1)
		
		# Compose a selection command which will select all atoms which are 
		sel_string = selection + " & ! b < " + str(lower)
		
		if(i < num_colors - 1):
			sel_string += " & b < " + str(upper)
		else:
			sel_string += " & ! b > " + str(upper)
		
		# Select the atoms
		cmd.select(group, sel_string) 

		# Create a new color
		color_name = "color_" + str(i+1)
		cmd.set_color(color_name, color)

		# color them
		cmd.color(color_name, group)


	# Create new color for insufficient sequences
	insuf_color = [0.75, 0.75, 0.58823529]
	cmd.set_color("insufficient_color", insuf_color)

	# color atoms with B-factor of 10 using the new color
	cmd.select("insufficient", selection + " & b = 0")
	cmd.color("insufficient_color", "insufficient")

# This is required to make command available in PyMOL 
cmd.extend("color_identity", color_identity)

color_identity()
