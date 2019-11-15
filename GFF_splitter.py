#!/usr/bin/env python3

import pandas as pd
import numpy as np
from os import listdir, remove
from os.path import isfile, join
import argparse
import sys

## USER ARGUMENT PARSING ##

## USER ARGUMENT PARSING ##
# take input from command line
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-g', '--GFF', type=str, help='GFF file')
args = parser.parse_args()
# check that a CG.output file has been given
if args.GFF:
	input_GFF = args.GFF
else:
	print("ERROR: no GFF file (-g) specified")
	sys.exit(0)

def GFF_splitter(GFF):
	"""
	A function that takes a GFF file and splits its contents by chromosome/contig into a set of new GFF files.
	It also writes a file chromosomes.txt which records the name of each chromosome in the GFF, one name per line
	"""
	## GFF SPLITTING ##
	# read in GFF
	features = pd.read_csv(GFF, sep = "\t", names = ["Chromosome", "Program", "Feature", "Start", "End", "Score", "Strand", "BLANK", "Metadata"])
	# sort by chromosome, start coordinate and feature name, and reset index
	features = features.sort_values(by = ["Chromosome", "Start", "Feature"], ascending = [True, True, False]).reset_index(drop = True)
	# output each chromosome's GFF entries into a separate GFF, with the filename as the name of that chromosome
	features_by_chromosome = features.groupby("Chromosome")
	for (chromosome, features_df) in features_by_chromosome:
		filename = chromosome + ".gff"
		features_df.to_csv(filename, sep = "\t", index = False, header = False)
	print("GFF file split into subfiles")
	
	## CHROMOSOME NAMES ##
	# find the chromosome names that this GFF file contained
	chromosome_names = sorted([chromosome for (chromosome, features_df) in features_by_chromosome])
	# output the chromosome names to file, one name per line, so that this can be used as input to GNU parallels 
	chromosomes_out = open("chromosomes.txt", "w")
	for i in chromosome_names:
		chromosomes_out.write(i + "\n")
	chromosomes_out.close()
	return()

GFF_splitter(GFF = input_GFF)
