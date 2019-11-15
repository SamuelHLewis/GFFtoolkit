import pandas as pd
import numpy as np
import argparse
import sys
import re

## USER ARGUMENT PARSING ##
# take input from command line
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-g', '--GFF', type=str, help='GFF file')
parser.add_argument('-l', '--sharedlabel', type=str, help='Label of field in GFF file that is shared between genes, exons and CDS')
args = parser.parse_args()
# check that a shared label has been given
if args.sharedlabel:
	input_sharedlabel = args.sharedlabel
else:
	print("ERROR: no shared label (--sharedlabel) specified")
	sys.exit(0)
# assign the GFF filename
if args.GFF:
	input_GFF = args.GFF
else:
	input_GFF = "No GFF"

def gene_exon_CDS_renamer(GFF, sharedlabel):
	"""
	A function that takes a GFF file, renames genes, exons and CDS to have the same name (as specified by the sharedlabel argument), and outputs a new GFF
	"""
	# read in GFF
	features = pd.read_csv(GFF, sep = "\t", names = ["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Metadata"])
	# remove any entries that were commented out
	features = features[features.Chromosome.str[0] != "#"]
	features.reset_index(inplace = True, drop = True)
#	print("GFF parsed")
#	print("After GFF reading:")
#	print(features.head())
	# define the search pattern for the shared name (which is identical between the CDS, exon and gene annotation files)
	shared_search_pattern = re.escape(sharedlabel) + r"=([\w\d_\-\:\.]+);"
	# extract the CDS, exon and gene annotations from the GFF
	CDS = pd.DataFrame(features[features["Feature"] == "CDS"])
	features.drop(CDS.index.tolist(), inplace = True)
	CDS.reset_index(inplace = True, drop = True)
	exons = pd.DataFrame(features[features["Feature"] == "exon"])
	features.drop(exons.index.tolist(), inplace = True)
	exons.reset_index(inplace = True, drop = True)
	genes = pd.DataFrame(features[features["Feature"] == "gene"])
	features.drop(genes.index.tolist(), inplace = True)
	genes.reset_index(inplace = True, drop = True)
	# create new columns in CDS, exon and gene dataframes for the shared name
	CDS["Shared_Name"] = CDS["Metadata"].str.extract(shared_search_pattern, flags=re.IGNORECASE, expand=False)
	exons["Shared_Name"] = exons["Metadata"].str.extract(shared_search_pattern, flags=re.IGNORECASE, expand=False)
	genes["Shared_Name"] = genes["Metadata"].str.extract(shared_search_pattern, flags=re.IGNORECASE, expand=False)
#	print("After shared name extraction")
#	print(CDS.head())
#	print(exons.head())
#	print(genes.head())
	# screen out any gene that doesn't have any CDS annotations (ie noncoding loci) from the gene and exon dataframes
	genes = genes[genes.Shared_Name.isin(CDS["Shared_Name"])]
	exons = exons[exons.Shared_Name.isin(CDS["Shared_Name"])]
#	print("After non-coding loci removal")
#	print(exons.head())
#	print(genes.head())
	# set the metadata for the genes, exons and CDS as the shared name, and drop the shared name column
	CDS["Metadata"] = CDS["Shared_Name"]
	CDS.drop(["Shared_Name"], axis = 1, inplace = True)
	exons["Metadata"] = exons["Shared_Name"]
	exons.drop(["Shared_Name"], axis = 1, inplace = True)
	genes["Metadata"] = genes["Shared_Name"]
	genes.drop(["Shared_Name"], axis = 1, inplace = True)
#	print("After metadata respecification")
#	print(CDS.head())
#	print(exons.head())
#	print(genes.head())
	# add the CDS, exon and gene annotations back into the GFF
	features = features.append(CDS, ignore_index = True)
	features = features.append(exons, ignore_index = True)
	features = features.append(genes, ignore_index = True)
	# respecify the start and end columns as integers
	features["Start"] = features["Start"].astype(int)
	features["End"] = features["End"].astype(int)
#	print("Final GFF file:")
#	print(features.head())
	features.sort_values(by = ["Chromosome", "Metadata", "Start"], inplace = True)
	features.to_csv("renamed.gff", header = False, index = False, sep = "\t")
	print("Renamed GFF has been written to renamed.gff")
	return()

gene_exon_CDS_renamer(GFF = input_GFF, sharedlabel = input_sharedlabel)

