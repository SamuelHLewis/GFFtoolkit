import os
import sys
import argparse
import re
import pandas as pd
import numpy as np

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-i', '--interpro', type=str, help='interpro output file (must be in tab-separated format)')
parser.add_argument('-g', '--gff', type=str, help='GFF file (the same GFF that was used as input for interpro)')
args = parser.parse_args()
# interpro file verification
interpro_file = args.interpro
if interpro_file is None:
	print('ERROR: no interpro output file (-i) specified')
	sys.exit(0)
# GFF file verification
GFF_file = args.gff
if GFF_file is None:
	print('ERROR: no GFF file (-g) specified')
	sys.exit(0)

def TE_relabeller(interpro, GFF):
	"""
	A function that takes a GFF file and an interpro output file, relabels TEs to give each a unique name, and outputs a GFF file with relabelled TE annotations.
	The unique names are of the form: Chromosome;Start;End;Strand;TE_name;TE_domain;TE_family;TE_class
	The TE domains are: PF03184 (DDE_1), PF02914 (DDE_2), PF13358 (DDE_3), PF03732 (Retrotrans_gag), PF00665 (rve) & PF00077 (RVP)
	"""
	# read in domain file
	domains = pd.read_csv(interpro, sep = "\t", names = ["interpro_name", "ignore1", "ignore2", "Domain", "Pfam_Domain", "ignore3", "ignore4", "ignore5", "ignore6", "ignore7", "ignore8", "ignore9", "ignore10", "ignore11", "ignore12", "ignore13", "ignore14", "ignore15", "ignore16", "ignore17"], dtype = "str")
	domains["Pfam_Domain"] = domains["Pfam_Domain"].str.replace(r"\.[0-9]*", "")
#	print("Domains parsed:")
#	print(domains.head())
	# read in feature file
	features = pd.read_csv(GFF, sep = "\t", names = ["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Metadata"], dtype = {"Chromosome":"str", "Program":"str", "Feature":"str", "Start":np.int32, "End":"str", "BLANK":"str", "Strand":"str", "Score":"str", "Metadata":"str"})
	# remove any entries that were commented out
	features = features[features.Chromosome.str[0] != "#"]
	features.reset_index(inplace = True, drop = True)
#	print("Features parsed:")
#	print(features.head())
	# create a new column for the start coordinate as assigned by gffread and used in interpro TE name
	features["Interpro_Start"] = features["Start"] - 1
	# drop all columns from domain file except the TE name and its domain
	domains = domains[["interpro_name", "Pfam_Domain"]]
#	print("Irrelevant data dropped from domains")
#	print(domains.head())
	# collapse multiple entries for the same TE into a single underscore-delimited entry
	domains_collapsed = pd.DataFrame(domains.groupby("interpro_name")["Pfam_Domain"].agg(lambda x:"_".join(x.unique())))
	domains_collapsed.reset_index(inplace = True)
#	print("Domains collapsed to unique entries with combined domain assignments")
#	print(domains_collapsed.head())
	# drop all rows that don't have a TE domain
	TE_domains = ["PF03184", "PF02914", "PF13358", "PF03732", "PF00665", "PF00077"]
	domains = domains_collapsed[domains_collapsed.Pfam_Domain.isin(TE_domains)]
#	print("Interpro entries for TE domains:")
#	print(domains.head())
	# adjust gene names in GFF to match those created by interpro
	features["interpro_name"] = features["Chromosome"] + ":" + features["Interpro_Start"].apply(str) + "-" + features["End"] + "(" + features["Strand"] + ")"
#	print(features.head())
	# add domain to each TE
	features = pd.merge(features, domains, how = "outer")
	print(features.shape[0] - features.isnull().sum().sum(), "TEs have at least 1 domain")
	# replace NaN with NA
	features.replace(np.nan, "NA", regex = True, inplace = True)
#	print("Domains merged with feature annotations")
	# change name of each TE to form Chromosome;Start;End;Strand;TE_name;TE_domain;TE_family;TE_class
	TE_new_names = []
	for i in features.itertuples():
		TE_new_names.append(";".join([i.Chromosome, str(i.Start), i.End, i.Strand, i.Metadata.split(";")[0], i.Pfam_Domain, i.Metadata.split(";")[1], i.Metadata.split(";")[2]]))
	features["Metadata"] = TE_new_names
	features.drop(["Interpro_Start", "interpro_name", "Pfam_Domain"], inplace = True, axis = 1)
#	print("TEs relabelled")
#	print(features.head())
#	# output to file
	features.to_csv("TE_relabelled.gff", sep = "\t", index = False, header = False)
	print("Relabelled TE GFF file written to TE_relabelled.gff")
	return

TE_relabeller(interpro = interpro_file, GFF = GFF_file)
