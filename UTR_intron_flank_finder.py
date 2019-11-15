#!/usr/bin/env python3

import pandas as pd
import numpy as np
from os import listdir, remove
from os.path import isfile, join
import argparse
import sys

## USER ARGUMENT PARSING ##
# take input from command line
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-g', '--GFF', type=str, help='GFF file')
args = parser.parse_args()
# check that a GFF file has been given
if args.GFF:
	input_GFF = args.GFF
else:
	print("ERROR: no GFF file (--GFF) specified")
	sys.exit(0)

def UTR_intron_flank_finder(GFF, flank_length = 1000):
	"""
	A function that takes a GFF file and adds entries for UTRs, introns in UTRs, and introns between CDS
	"""
	## GFF IMPORT AND PREPROCESSING ##
	# read in GFF
	features = pd.read_csv(GFF, sep = "\t", names = ["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Metadata"])
#	print("GFF read")
	# make a new dataframe holding the coding span for each gene
	coding_spans = pd.DataFrame(features[features["Feature"] == "CDS"].groupby("Metadata")["Start"].min()).join(pd.DataFrame(features[features["Feature"] == "CDS"].groupby("Metadata")["End"].max())).join(pd.DataFrame(features[features["Feature"] == "CDS"].groupby("Metadata")["Strand"].first()))
#	print("Coding spans calculated")
#	print(coding_spans.head())

	## FIND UTRs ##
	total_features = features.shape[0]
	for feature in features.itertuples():
		# for each exon that sits in a gene with coding regions
		if feature.Feature == "exon" and feature.Metadata in coding_spans.index.tolist():
			# look up coding span start and end for this gene
			coding_start = coding_spans.loc[feature.Metadata]["Start"]
			coding_end = coding_spans.loc[feature.Metadata]["End"]
			# deal with lefthand UTRs that aren't contiguous with CDS
			if feature.Start < coding_start and feature.End < coding_start:
				if feature.Strand == "+":
					features.iat[feature.Index, 2] = "five_prime_UTR"
				elif feature.Strand == "-":
					features.iat[feature.Index, 2] = "three_prime_UTR"
			# deal with lefthand UTRs that are contiguous with CDS
			elif feature.Start < coding_start and feature.End > coding_start:
				features.iat[feature.Index, 4] = coding_start - 1
				if feature.Strand == "+":
					features.iat[feature.Index, 2] = "five_prime_UTR"
				elif feature.Strand == "-":
					features.iat[feature.Index, 2] = "three_prime_UTR"
			# deal with righthand UTRs that aren't contiguous with CDS
			elif feature.Start > coding_end and feature.End > coding_end:
				if feature.Strand == "+":
					features.iat[feature.Index, 2] = "three_prime_UTR"
				elif feature.Strand == "-":
					features.iat[feature.Index, 2] = "five_prime_UTR"
			# deal with righthand UTRs that are contiguous with CDS
			elif feature.Start < coding_end and feature.End > coding_end:
				features.iat[feature.Index, 3] = coding_end + 1
				if feature.Strand == "+":
					features.iat[feature.Index, 2] = "three_prime_UTR"
				elif feature.Strand == "-":
					features.iat[feature.Index, 2] = "five_prime_UTR"
			# change the Program field for all former exons to reflect that their new content was produced by this script
			features.iat[feature.Index, 1] = "UTR_intron_flank_finder.py"
	# the only exons left now will overlap entirely with CDS or be exons from non-coding features, so remove them and resort the gff
	features = features[features["Feature"] != "exon"].sort_values(by = ["Chromosome", "Start"]).reset_index(drop = True)
#	print("UTRs found")
#	print(features)
	## FIND INTRONS ##
	# first, find introns between CDS
	for gene in coding_spans.itertuples():
		# grab the CDS entries for this gene
		CDS = features[(features["Feature"] == "CDS") & (features["Metadata"] == gene.Index)]
		# for each CDS except the last, define an intron entry based on the end of the CDS +1 and the start of the next CDS -1
		CDS_data = list(CDS.itertuples())
		for i in range(len(CDS_data)-1):
			intron_start = CDS_data[i].End + 1
			intron_end = CDS_data[i+1].Start - 1
			intron_entry = pd.DataFrame({"Chromosome":CDS_data[i].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron", "Start":intron_start, "End":intron_end, "BLANK":".", "Strand":CDS_data[i].Strand, "Score":".", "Metadata":CDS_data[i].Metadata}, index = [0])
			features = features.append(intron_entry, ignore_index = True)	
#	print("Introns between CDS found")
#	print(features.head())
	# next, find introns in UTRs
	# for each gene
	for CDS_span in coding_spans.itertuples():
		# calculate the left and right non-coding spans
		left_noncoding_start = features[(features["Feature"] == "gene") & (features["Metadata"] == CDS_span.Index)]["Start"].values[0]
		left_noncoding_end = CDS_span.Start - 1
#		print("The left non-coding region for", CDS_span.Index, "runs from", left_noncoding_start, "to", left_noncoding_end)
		right_noncoding_start = CDS_span.End + 1
		right_noncoding_end = features[(features["Feature"] == "gene") & (features["Metadata"] == CDS_span.Index)]["End"].values[0]
#		print("The right non-coding region for", CDS_span.Index, "runs from", right_noncoding_start, "to", right_noncoding_end)
		# grab the five-prime and three-prime UTR entries for this gene
		five_prime_UTRs = features[(features["Feature"] == "five_prime_UTR") & (features["Metadata"] == CDS_span.Index)]
#		print(five_prime_UTRs)
		three_prime_UTRs = features[(features["Feature"] == "three_prime_UTR") & (features["Metadata"] == CDS_span.Index)]
#		print(three_prime_UTRs)
		# for genes on the + strand that have at least 1 5'UTR
		if CDS_span.Strand == "+" and five_prime_UTRs.shape[0] > 0:
			# find introns in the left non-coding region
			# check if there is an intron between the start of the left non-coding region and the start of the leftmost five-prime UTR
			if left_noncoding_start < five_prime_UTRs.iloc[0]["Start"]:
				intron_entry = pd.DataFrame({"Chromosome":five_prime_UTRs.iloc[0].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_five_prime_UTR", "Start":left_noncoding_start, "End":five_prime_UTRs.iloc[0]["Start"]-1, "BLANK":".", "Strand":five_prime_UTRs.iloc[0].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
				features = features.append(intron_entry, ignore_index = True)
			# if there is more than one five-prime UTR, find introns between UTRs
			if five_prime_UTRs.shape[0] > 1:
				UTR_data = list(five_prime_UTRs.itertuples())
				for i in range(len(UTR_data)-1):
					intron_start = UTR_data[i].End + 1
					intron_end = UTR_data[i+1].Start - 1
					intron_entry = pd.DataFrame({"Chromosome":five_prime_UTRs.iloc[0].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_five_prime_UTR", "Start":intron_start, "End":intron_end, "BLANK":".", "Strand":five_prime_UTRs.iloc[0].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
					features = features.append(intron_entry, ignore_index = True)
			# check if there is an intron between the end of the rightmost five-prime UTR and the end of the left non-coding region
			if five_prime_UTRs.iloc[-1]["End"] < left_noncoding_end:
				intron_entry = pd.DataFrame({"Chromosome":five_prime_UTRs.iloc[-1].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_five_prime_UTR", "Start":five_prime_UTRs.iloc[-1]["End"]+1, "End":left_noncoding_end, "BLANK":".", "Strand":five_prime_UTRs.iloc[-1].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
				features = features.append(intron_entry, ignore_index = True)
		# for genes on the + strand that have at least 1 3'UTR
		if CDS_span.Strand == "+" and three_prime_UTRs.shape[0] > 0:
			# find introns in the right non-coding region
			# check if there is an intron between the start of the right non-coding region and the start of the leftmost three-prime UTR
			if right_noncoding_start < three_prime_UTRs.iloc[0]["Start"]:
				intron_entry = pd.DataFrame({"Chromosome":three_prime_UTRs.iloc[0].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_three_prime_UTR", "Start":right_noncoding_start, "End":three_prime_UTRs.iloc[0]["Start"]-1, "BLANK":".", "Strand":three_prime_UTRs.iloc[0].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
			# if there is more than one three-prime UTR, find introns between UTRs
			if three_prime_UTRs.shape[0] > 1:
				UTR_data = list(three_prime_UTRs.itertuples())
				for i in range(len(UTR_data)-1):
					intron_start = UTR_data[i].End + 1
					intron_end = UTR_data[i+1].Start - 1
					intron_entry = pd.DataFrame({"Chromosome":three_prime_UTRs.iloc[0].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_three_prime_UTR", "Start":intron_start, "End":intron_end, "BLANK":".", "Strand":three_prime_UTRs.iloc[0].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
					features = features.append(intron_entry, ignore_index = True)
			# check if there is an intron between the end of the rightmost three-prime UTR and the end of the right non-coding region
			if three_prime_UTRs.iloc[-1]["End"] < right_noncoding_end:
				intron_entry = pd.DataFrame({"Chromosome":three_prime_UTRs.iloc[-1].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_three_prime_UTR", "Start":three_prime_UTRs.iloc[-1]["End"]+1, "End":right_noncoding_end, "BLANK":".", "Strand":three_prime_UTRs.iloc[-1].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
				features = features.append(intron_entry, ignore_index = True)
		# for genes on the - strand with at least one 3'UTR
		if CDS_span.Strand == "-" and three_prime_UTRs.shape[0] > 0:
			# find introns in the left non-coding region
			# check if there is an intron between the start of the left non-coding region and the start of the leftmost three-prime UTR
			if right_noncoding_start < three_prime_UTRs.iloc[0]["Start"]:
				intron_entry = pd.DataFrame({"Chromosome":three_prime_UTRs.iloc[0].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_three_prime_UTR", "Start":left_noncoding_start, "End":three_prime_UTRs.iloc[0]["Start"]-1, "BLANK":".", "Strand":three_prime_UTRs.iloc[0].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
				features = features.append(intron_entry, ignore_index = True)
			# if there is more than one three-prime UTR, find introns between UTRs
			if three_prime_UTRs.shape[0] > 1:
				UTR_data = list(three_prime_UTRs.itertuples())
				for i in range(len(UTR_data)-1):
					intron_start = UTR_data[i].End + 1
					intron_end = UTR_data[i+1].Start - 1
					intron_entry = pd.DataFrame({"Chromosome":three_prime_UTRs.iloc[0].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_three_prime_UTR", "Start":intron_start, "End":intron_end, "BLANK":".", "Strand":three_prime_UTRs.iloc[0].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
					features = features.append(intron_entry, ignore_index = True)
			# check if there is an intron between the end of the rightmost three-prime UTR and the end of the left non-coding region
			if three_prime_UTRs.iloc[-1]["End"] < left_noncoding_end:
				intron_entry = pd.DataFrame({"Chromosome":three_prime_UTRs.iloc[-1].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_three_prime_UTR", "Start":three_prime_UTRs.iloc[-1]["End"]+1, "End":left_noncoding_end, "BLANK":".", "Strand":three_prime_UTRs.iloc[-1].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
				features = features.append(intron_entry, ignore_index = True)
		# for genes on the - strand that have at least 1 3'UTR
		if CDS_span.Strand == "-" and five_prime_UTRs.shape[0] > 0:
			# find introns in the right non-coding region
			# check if there is an intron between the start of the right non-coding region and the start of the leftmost five-prime UTR
			if right_noncoding_start < five_prime_UTRs.iloc[0]["Start"]:
				intron_entry = pd.DataFrame({"Chromosome":five_prime_UTRs.iloc[0].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_five_prime_UTR", "Start":right_noncoding_start, "End":five_prime_UTRs.iloc[0]["Start"]-1, "BLANK":".", "Strand":five_prime_UTRs.iloc[0].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
			# if there is more than one five-prime UTR, find introns between UTRs
			if five_prime_UTRs.shape[0] > 1:
				UTR_data = list(five_prime_UTRs.itertuples())
				for i in range(len(UTR_data)-1):
					intron_start = UTR_data[i].End + 1
					intron_end = UTR_data[i+1].Start - 1
					intron_entry = pd.DataFrame({"Chromosome":five_prime_UTRs.iloc[0].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_five_prime_UTR", "Start":intron_start, "End":intron_end, "BLANK":".", "Strand":five_prime_UTRs.iloc[0].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
					features = features.append(intron_entry, ignore_index = True)
			# check if there is an intron between the end of the rightmost five-prime UTR and the end of the right non-coding region
			if five_prime_UTRs.iloc[-1]["End"] < right_noncoding_end:
				intron_entry = pd.DataFrame({"Chromosome":five_prime_UTRs.iloc[-1].Chromosome, "Program":"UTR_intron_flank_finder.py", "Feature":"intron_five_prime_UTR", "Start":five_prime_UTRs.iloc[-1]["End"]+1, "End":right_noncoding_end, "BLANK":".", "Strand":five_prime_UTRs.iloc[-1].Strand, "Score":".", "Metadata":CDS_span.Index}, index = [0])
				features = features.append(intron_entry, ignore_index = True)
#	print("Introns in UTRs found")
	# sort the gff after adding introns
	features = features.sort_values(by = ["Chromosome", "Start"]).reset_index(drop = True)
#	print("After adding introns:")
#	print(features.head())
		
	## FLANKING POSITION COMPUTATION ##
	# filter out anything that isn't a gene, TE, rRNA or tRNA annotation
	features_for_flanks = features[(features["Feature"] == "gene") | (features["Feature"] == "TE") | (features["Feature"] == "rRNA") | (features["Feature"] == "tRNA")]
	# add a new entry to the features dataframe for each upstream and downstream flank
	for feature in features_for_flanks.itertuples():
		if feature.Strand == "+":
			# only add upstream flanks where they don't overhang the start of the chromosome
			if feature.Start > flank_length:
				upstream_flank = pd.DataFrame({"Chromosome": feature.Chromosome, "Program": "UTR_intron_flank_finder.py", "Feature": feature.Feature + "_upstream_flank", "Start": feature.Start - flank_length, "End": feature.Start - 1, "BLANK": ".", "Strand": "+", "Score": ".", "Metadata": feature.Metadata}, index = [0])
				features = features.append(upstream_flank, ignore_index = True)
#				print("Upstream flank for + strand feature appended")
			downstream_flank = pd.DataFrame({"Chromosome": feature.Chromosome, "Program": "UTR_intron_flank_finder.py", "Feature": feature.Feature + "_downstream_flank", "Start": feature.End + 1, "End": feature.End + flank_length, "BLANK": ".", "Strand": "+", "Score": ".", "Metadata": feature.Metadata}, index = [0])
			features = features.append(downstream_flank, ignore_index = True)
#			print("Downstream flank for + strand feature appended")
		elif feature.Strand == "-":
			upstream_flank = pd.DataFrame({"Chromosome": feature.Chromosome, "Program": "UTR_intron_flank_finder.py", "Feature": feature.Feature + "_upstream_flank", "Start": feature.End + 1, "End": feature.End + flank_length, "BLANK": ".", "Strand": "-", "Score": ".", "Metadata": feature.Metadata}, index = [0])
			features = features.append(upstream_flank, ignore_index = True)
#			print("Upstream flank for - strand feature appended")
			# only add downstream flanks where they don't overhang the start of the chromosome
			if feature.Start > flank_length:
				downstream_flank = pd.DataFrame({"Chromosome": feature.Chromosome, "Program": "UTR_intron_flank_finder.py", "Feature": feature.Feature + "_downstream_flank", "Start": feature.Start - flank_length, "End": feature.Start - 1, "BLANK": ".", "Strand": "-", "Score": ".", "Metadata": feature.Metadata}, index = [0])
				features = features.append(downstream_flank, ignore_index = True)
#				print("Downstream flank for - strand feature appended")
	# sort the gff after adding flanking regions
	features = features.sort_values(by = ["Chromosome", "Start"]).reset_index(drop = True)
	features = features[["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Metadata"]]
#	print("Flanking regions found")
#	print(features.head())

	## FILE OUTPUT ##
	filename = input_GFF.replace(".gff", ".with_introns_UTRs_flanks.gff")
	features.to_csv(filename, sep = "\t", index = False, header = False)
	print("GFF file with introns, UTRs and flanks written to", filename)
	return()

UTR_intron_flank_finder(GFF = input_GFF)

