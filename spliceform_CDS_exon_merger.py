import pandas as pd
import numpy as np
import argparse
import sys

## USER ARGUMENT PARSING ##
# take input from command line
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-g', '--GFF', type=str, help='GFF file')
args = parser.parse_args()
# assign the GFF filename
if args.GFF:
	input_GFF = args.GFF
else:
	input_GFF = "No GFF"

def spliceform_CDS_exon_merger(GFF_df):
	"""
	A function that takes a DataFrame parsed from a GFF file, merges any overlapping CDS annotations for each gene into a single CDS annotation, and outputs a tuple of the DataFrame (usually merged, but not always) and the status of whether any merging occured (boolean)
	"""
	features = GFF_df
	# separate out CDS annotations and sort by gene name
	CDS = features[features["Feature"] == "CDS"]
	features.drop(CDS.index.tolist(), inplace = True)
	CDS.sort_values(["Metadata", "Start"], inplace = True)
	CDS = CDS.reset_index(drop = True)
#	print(CDS)
	# separate out exon annotations and sort by gene name
	exons = features[features["Feature"] == "exon"]
	features.drop(exons.index.tolist(), inplace = True)
	exons.sort_values(["Metadata", "Start"], inplace = True)
	exons = exons.reset_index(drop = True)
#	print(exons)
	# find overlapping CDS annotions
	CDS_to_add = pd.DataFrame(columns = ["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Metadata"])
	CDS_to_delete = []
	for feature in CDS.itertuples():
		if feature.Index == 0:
			continue
#		print("Index of this feature =", feature.Index)
#		print("Start of this feature =", feature.Start)
#		print("End of feature above this one:")
#		print(CDS.iloc[feature.Index - 1]["End"])
		if int(feature.Start) <= int(CDS.iloc[int(feature.Index) - 1].End) and feature.Metadata == CDS.iloc[int(feature.Index) - 1].Metadata:
#			print("Found overlapping CDS:")
#			print(CDS.iloc[feature.Index - 1])
#			print(CDS.iloc[feature.Index])
			new_start = CDS.iloc[int(feature.Index) - 1].Start
			# check which end of the overlapping feature is the biggest, and set this as the new end position
			if feature.End >= CDS.iloc[int(feature.Index) - 1].End:
				new_end = feature.End
			else:
				new_end = CDS.iloc[int(feature.Index) - 1].End
			merged_entry = pd.DataFrame({
				"Chromosome": feature.Chromosome,
				"Program": feature.Program,
				"Feature": feature.Feature,
				"Start": str(new_start),
				"End": str(new_end),
				"BLANK": ".",
				"Strand": feature.Strand,
				"Score": str(feature.Score),
				"Metadata": feature.Metadata
			}, index = [0])
#			print("New merged entry:")
#			print(merged_entry)
			CDS_to_add = CDS_to_add.append(merged_entry)
			CDS_to_delete.extend([feature.Index - 1, feature.Index])
	if len(CDS_to_delete) > 0:
		# write names of overlapping CDS to file
		CDS.iloc[CDS_to_delete].to_csv("overlapping_CDS.txt", header = False, index = False, sep = "\t", mode = "a")
#		print("Overlapping CDS written to overlapping_CDS.txt")
		# delete overlapping CDS annotations
		CDS.drop(CDS_to_delete, inplace = True)
		# add merged CDS annotations
		CDS = CDS.append(CDS_to_add)
		CDS.sort_values(by = ["Chromosome", "Start"], inplace = True)
		CDS.reset_index(drop = True, inplace = True)
#		print(len(CDS_to_delete), "overlapping CDS merged")
		# set the merged CDS flag to True
		merged_CDS_flag = True
	else:
#		print("No overlapping CDS detected")
		# set the merged CDS flag to False
		merged_CDS_flag = False
#	print("After deleting overlapping CDS and adding merged CDS:")
#	print(CDS)
	# find overlapping exon annotions
	exons_to_add = pd.DataFrame(columns = ["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Metadata"])
	exons_to_delete = []
	for feature in exons.itertuples():
		if feature.Index == 0:
			continue
#		print("Index of this feature =", feature.Index)
#		print("Start of this feature =", feature.Start)
#		print("End of feature above this one:")
#		print(exons.iloc[feature.Index - 1]["End"])
		if int(feature.Start) <= int(exons.iloc[int(feature.Index) - 1].End) and feature.Metadata == exons.iloc[int(feature.Index) - 1].Metadata:
#			print("Found overlapping exon:")
#			print(exons.iloc[feature.Index - 1])
#			print(exons.iloc[feature.Index])
			new_start = exons.iloc[int(feature.Index) - 1].Start
			# check which end of the overlapping feature is the biggest, and set this as the new end position
			if feature.End >= exons.iloc[int(feature.Index) - 1].End:
				new_end = feature.End
			else:
				new_end = exons.iloc[int(feature.Index) - 1].End
			merged_entry = pd.DataFrame({
				"Chromosome": feature.Chromosome,
				"Program": feature.Program,
				"Feature": feature.Feature,
				"Start": str(new_start),
				"End": str(new_end),
				"BLANK": ".",
				"Strand": feature.Strand,
				"Score": str(feature.Score),
				"Metadata": feature.Metadata
			}, index = [0])
#			print("New merged entry:")
#			print(merged_entry)
			exons_to_add = exons_to_add.append(merged_entry)
			exons_to_delete.extend([feature.Index - 1, feature.Index])
	if len(exons_to_delete) > 0:
		# write names of overlapping exons to file
		exons.iloc[exons_to_delete].to_csv("overlapping_exons.txt", header = False, index = False, sep = "\t", mode = "a")
#		print("Overlapping exons written to overlapping_exons.txt")
		# delete overlapping exon annotations
		exons.drop(exons_to_delete, inplace = True)
		# add merged exon annotations
		exons = exons.append(exons_to_add)
		exons.sort_values(by = ["Chromosome", "Start"], inplace = True)
		exons.reset_index(drop = True, inplace = True)
#		print(len(exons_to_delete), "overlapping exons merged")
		# set the merged exons flag as True
		merged_exons_flag = True
	else:
#		print("No overlapping exons detected")
		merged_exons_flag = False
#	print("After deleting overlapping exons and adding merged exons:")
#	print(exons)
	# add CDS and exon annotations back into all features df
	features = features.append(CDS)
	features = features.append(exons)
	# respecify the start and end columns as integers
	features["Start"] = features["Start"].astype(int)
	features["End"] = features["End"].astype(int)
	# sort the features by chromosome name, gene name, start coordinate, and feature name
	features = features.sort_values(by = ["Chromosome", "Metadata", "Start", "Feature"], ascending = [True, True, True, False])
	features.reset_index(drop = True, inplace = True)
	# reorder columns to reflect GFF format
	features = features[["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Metadata"]]
#	print("Final GFF file:")
#	print(features)
#	features.to_csv("merged.gff", header = False, index = False, sep = "\t")
#	print("Merged GFF written to merged.gff")
	# set the overall merge flag as True or False if either CDS or exons have been merged
	if merged_CDS_flag == True or merged_exons_flag == True:
		return((features, True))
	else:
		return((features, False))

# read in GFF
features = pd.read_csv(input_GFF, sep = "\t", names = ["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Metadata"])
# remove any entries that were commented out
features = features[features.Chromosome.str[0] != "#"]
features.reset_index(inplace = True, drop = True)
print("GFF parsed")
# print("After GFF reading:")
# print(features)

# write blank files for gene names of merged CDS and exons (so that names can be appended to these after each merge cycle)
with open("overlapping_CDS.txt", "wt") as overlapping_CDS_outfile:
	overlapping_CDS_outfile.write("Cycle 1:\n")
with open("overlapping_exons.txt", "wt") as overlapping_exons_outfile:
	overlapping_exons_outfile.write("Cycle 1:\n")

# do the first merge
merged_GFF = spliceform_CDS_exon_merger(GFF_df = features)
merge_count = 1
more_merges_needed = merged_GFF[1]
print("Merge cycle 1 complete")

# do further merges until all CDS/exons are merged
while more_merges_needed == True:
	merge_count += 1
	with open("overlapping_CDS.txt", "a") as overlapping_CDS_outfile:
		overlapping_CDS_outfile.write("Cycle " + str(merge_count) + ":\n")
	with open("overlapping_exons.txt", "a") as overlapping_exons_outfile:
		overlapping_exons_outfile.write("Cycle " + str(merge_count) + ":\n")
	merged_GFF = spliceform_CDS_exon_merger(GFF_df = merged_GFF[0])
	more_merges_needed = merged_GFF[1]
	print("Merge cycle", merge_count, "complete")

# output merged GFF to file
merged_GFF[0].to_csv("merged.gff", header = False, index = False, sep = "\t")
print("No remaining merges needed after merge cycle", merge_count - 1)
print("Merged GFF written to merged.gff")

