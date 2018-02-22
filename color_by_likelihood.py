#!/usr/bin/env python
#
# color_by_likelihood.py v1 2018-02-01

'''color_by_likelihood.py v1.1 2018-02-22

in PyMOL, use by:

run color_by_likelihood.py
'''

def color_likelihood(selection="all", whitebg=False):
	'''function to color atoms/residues by the temperatureFactor field'''

	# each color is RGB triplet, values ranging from 0 to 1, where black is 0 and white is 1
	# colors should be:
	# for values [ -1 gray
	#            0 green   1   2
	#            3 sea-green   4   5 teal
	#            6 verdigris 7 mud
	#            8 clay   9 pink   a
	#            b rose   c        d
	#            e red    f
	#            16 blue ]

	colors = [ [0.39, 0.39, 0.39] ,
             [0.11,0.57,0.06] , [0.12,0.63,0.21] , [0.19,0.74,0.34] ,
             [0.18,0.77,0.48] , [0.16,0.82,0.67] , [0.26,0.77,0.64] ,
             [0.39,0.72,0.60] , [0.42,0.49,0.42] ,
             [0.49,0.36,0.40] , [0.72,0.33,0.56] , [0.76,0.26,0.56] ,
             [0.78,0.24,0.48] , [0.77,0.22,0.44] , [0.77,0.20,0.33] ,
             [0.77,0.18,0.22] , [0.77,0.16,0.14] ,
             [0.12,0.44,0.88] ]

	if whitebg: # replaces colors -1, 7, 8
		colors = [ [0.95,1.00,0.50] ,
             [0.11,0.57,0.06] , [0.12,0.63,0.21] , [0.19,0.74,0.34] ,
             [0.18,0.77,0.48] , [0.16,0.82,0.67] , [0.26,0.77,0.64] ,
             [0.39,0.72,0.60] , [0.82,0.95,0.87] ,
             [0.93,0.80,0.78] , [0.72,0.33,0.56] , [0.76,0.26,0.56] ,
             [0.78,0.24,0.48] , [0.77,0.22,0.44] , [0.77,0.20,0.33] ,
             [0.77,0.18,0.22] , [0.77,0.16,0.14] ,
             [0.12,0.23,0.66] ]

	# v1 values of blue to teal
    #         [0.00,0.22,0.72] , [0.05,0.38,0.71] , [0.07,0.47,0.70] ,
    #         [0.11,0.63,0.69] , [0.16,0.82,0.67] , [0.26,0.77,0.64] ,
	# v1 constant color as green
    #         [0.21,0.66,0.00] ]

	# colors for:
	# 7d7d00 dark yellow    b28b71 sand    640064 burgundy
	# 966400 orange    785583 lilac    4b6419 brown
	gapcolors = [ [0.49,0.49,0.00] , [0.70,0.55,0.44] , [0.39,0.00,0.39] ,
                 [0.59,0.39,0.00] , [0.47,0.33,0.51] , [0.29,0.39,0.10] ] 

	groupnames = [ "gaps", 
                  "t2max", "-5", "-4", "-3", "-2", "-1", "-0.5", "-0.1",
                  "+0.1", "+0.5", "+1", "+2", "+3", "+4", "+5", "t1max",
                   "const" ]

	bin_size = 1

	# iterate through temperatureFactor value in PyMOL
	# atoms with each value are selected and colored with the following steps
	for i in range(18):

		lower = -1 + i
		upper = lower + 0.99

		# Print out B-factor limits and the color for this group
		print lower, " - ", upper, " = ", colors[i]

		# Define a unique name for the atoms which fall into this group
		groupname = groupnames[i] + "_group_s" + str(i+1)

		# Compose a selection command which will select all atoms with some beta value
		sel_string = selection + " & ! b < " + str(lower)

		if(i < 18 - 1):
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
	cmd.select("insufficient", selection + " & b > 16")
	cmd.color("insufficient_color", "insufficient")

	print "to use white background color scheme, type 'color_likelihood(whitebg=True)'"

# This is required to make command available in PyMOL 
cmd.extend("color_likelihood", color_likelihood)

color_likelihood()
