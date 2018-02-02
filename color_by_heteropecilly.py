#!/usr/bin/env python
#
# color_by_heteropecilly.py v1 2017-10-06

'''color_by_heteropecilly.py v1.1 2018-02-02

in PyMOL, use by:

run color_by_heteropecilly.py
'''

def color_heteropecilly(selection="all", num_colors=13, MINIMUM=-1.0, MAXIMUM=11.0):
	'''function to color atoms/residues by the temperatureFactor field'''

	# each color is RGB triplet, values ranging from 0 to 1, where black is 0 and white is 1
	# colors should be:
	# gray                   for no data or gap
    # blue-gray, cobalt blue, blue,
    # light blue, pale blue, pale pink,
    # pink, medium red, more red,
    # red                    for most heteropecillious
    # light green, green     for mostly constant (holozoa missing) and constant

	colors = [ [0.39, 0.39, 0.39] , 
              [0.38,0.41,0.55] , [0.37,0.43,0.73] , [0.36,0.45,0.85] ,
              [0.51,0.56,0.86] , [0.71,0.71,0.86] , [0.86,0.71,0.71] ,
              [0.87,0.50,0.50] , [0.87,0.36,0.36] , [0.87,0.22,0.22] , 
              [0.87,0.00,0.00] ,
              [0.44,0.78,0.45] , [0.21,0.66,0.00] ]
	groupnames = ["gaps", 
                "hp_group_0-10", "hp_group_10-20", "hp_group_20-30",
                "hp_group_30-40", "hp_group_40-50", "hp_group_50-60",
                "hp_group_60-70", "hp_group_70-80", "hp_group_80-90",
                "hp_group_90-99", "semi-const", "constant" ]

	bin_size = ((MAXIMUM - MINIMUM) + 1) / num_colors

	# iterate through temperatureFactor value in PyMOL
	# atoms with each value are selected and colored with the following steps
	for i in range(num_colors):

		lower = MINIMUM + i * bin_size
		upper = lower + bin_size - 0.01

		# Print out B-factor limits and the color for this group
		print lower, " - ", upper, " = ", colors[i]

		# Define a unique name for the atoms which fall into this group
		groupname = groupnames[i] + "_s" + str(i+1)

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
cmd.extend("color_heteropecilly", color_heteropecilly)

color_heteropecilly()
