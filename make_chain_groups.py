#!/usr/bin/env python
#
# make_chain_groups.py v1 2019-02-14

'''make_chain_groups.py  last modified 2019-02-14
    generate a PyMOL script to make selection groups for each chain

make_chain_groups.py -p 5ara.pdb > 5ara_chains.pml

    run script in the PyMOL console as:
@5ara_chains.pml
'''

import sys
import argparse

# DBREF specs from:
# http://www.wwpdb.org/documentation/file-format-content/format33/sect3.html#DBREF
#COLUMNS       DATA TYPE     FIELD              DEFINITION
#-----------------------------------------------------------------------------------
# 1 -  6       Record name   "DBREF "
# 8 - 11       IDcode        idCode             ID code of this entry.
#13            Character     chainID            Chain  identifier.
#15 - 18       Integer       seqBegin           Initial sequence number of the
#                                               PDB sequence segment.
#19            AChar         insertBegin        Initial  insertion code of the
#                                               PDB  sequence segment.
#21 - 24       Integer       seqEnd             Ending sequence number of the
#                                               PDB  sequence segment.
#25            AChar         insertEnd          Ending insertion code of the
#                                               PDB  sequence segment.
#27 - 32       LString       database           Sequence database name.
#34 - 41       LString       dbAccession        Sequence database accession code.
#43 - 54       LString       dbIdCode           Sequence  database identification code.
#56 - 60       Integer       dbseqBegin         Initial sequence number of the
#                                               database seqment.
#61            AChar         idbnsBeg           Insertion code of initial residue of the
#                                               segment, if PDB is the reference.
#63 - 67       Integer       dbseqEnd           Ending sequence number of the
#                                               database segment.
#68            AChar         dbinsEnd           Insertion code of the ending residue of
#                                               the segment, if PDB is the reference.

def make_chain_select_commands(pdbfile, wayout):
	'''read PDB file and return a dict where key is chain and value is sequence ID'''
	keepchains = {} # dict where key is chain and value is seqid
	chaintracker = {} # key is chains from ATOM records, value is True
	print >> sys.stderr, "# Reading chain info from PDB {}".format(pdbfile)
	for line in open(pdbfile,'r'):
		record = line[0:6].strip()
		# get relevant chains that match the sequence, in case of hetero multimers
		if record=="DBREF":
			proteinid = line[42:56].strip()
			chaintarget = line[12]
			keepchains[chaintarget] = proteinid
		elif record=="ATOM":
			chaintarget = line[21]
			chaintracker[chaintarget] = True
	if keepchains:
		for chain in sorted(keepchains.keys()):
			print >> wayout, "select {}__{}, chain {}".format( chain, keepchains[chain], chain )
		print >> sys.stderr, "# wrote 'select' commands for {} chains".format( len(keepchains) )
	else:
		print >> sys.stderr, "WARNING: NO DBREF RECORDS FOUND IN {}".format(pdbfile)
		if chaintracker:
			print >> sys.stderr, "ATOMS WERE FOUND FOR THE FOLLOWING CHAINS:\n{}".format( ",".join( sorted(chaintracker.keys() ) ) )
	# no return

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	args = parser.parse_args(argv)

	# make PyMOL script with color commands
	make_chain_select_commands(args.pdb, wayout)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
