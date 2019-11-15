#!/usr/bin/env python3

import os
import sys
import argparse
import re

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-i', '--input', type=str, help='input file (i.e. output file from tRNAscan)')
args = parser.parse_args()
# input file verification
tRNAscan_file = args.input
if tRNAscan_file is None:
	print('ERROR: no tRNAscan file specified')
	sys.exit(0)

##################################
## tRNAscan output file parsing ##
##################################
tRNA_chromosomes = []
tRNA_starts = []
tRNA_ends = []
tRNA_strands = []
tRNA_types = []
tRNA_anticodons = []
tRNA_scores = []
for line in open(tRNAscan_file, "r"):
	# skip headers
	if not line.startswith("Sequence") and not line.startswith("Name") and not line.startswith("--------"):
		# split lines on whitespace
		temp = line.split()
		# read start and end, and assign strand
		start = temp[2]
		end = temp[3]
		if int(start)<int(end):
			tRNA_start = start
			tRNA_end = end
			tRNA_strands.append("+")
		else:
			tRNA_start = end
			tRNA_end = start
			tRNA_strands.append("-")
		# split information for each line to lists
		tRNA_chromosomes.append(temp[0])
		tRNA_starts.append(tRNA_start)
		tRNA_ends.append(tRNA_end)
		tRNA_types.append(temp[4])
		tRNA_anticodons.append(temp[5])
		tRNA_scores.append(temp[8].strip("\n"))

###########################
## GFF output formatting ##
###########################
GFFoutput = ""
for i in range(len(tRNA_chromosomes)):
	GFFoutput += tRNA_chromosomes[i] + "\ttRNAscan\ttRNA\t" + tRNA_starts[i] + "\t" + tRNA_ends[i] + "\t" + tRNA_scores[i] + "\t" + tRNA_strands[i] + "\t.\tType=" + tRNA_types[i]  + ";Anticodon=" + tRNA_anticodons[i] + "\n"
# output the metadata for each name to a file
outfile=open("tRNAscan.gff","wt")
outfile.write(GFFoutput)
outfile.close()
print("GFF format output written to tRNAscan.gff")

