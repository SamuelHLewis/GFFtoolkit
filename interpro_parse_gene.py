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
parser.add_argument('-g', '--gff', type=str, help='GFF file (the same GFF that was used to generate proteins for interpro)')
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

def interpro_parse_gene(interpro, GFF):
	"""
	A function to rename GFF gene entries that have TE domains, as identified by interproscan.
	The TE domains are: PF03184 (DDE_1), PF02914 (DDE_2), PF13358 (DDE_3), PF03732 (Retrotrans_gag), PF00665 (rve) & PF00077 (RVP)
	"""
	# read in input files
	domains = pd.read_csv(interpro, sep = "\t", names = ["Gene", "MD5", "Length", "Algorithm", "Domain", "Description", "Start", "End", "Score", "Status", "Date", "return", "Interpro_Accession", "Interpro_Description"], dtype = "str")
	features = pd.read_csv(GFF, sep = "\t", names = ["Chromosome", "Program", "Feature", "Start", "End", "BLANK", "Strand", "Score", "Gene"], dtype = "str")
	# remove any entries that were commented out
	features = features[features.Chromosome.str[0] != "#"]
	features.reset_index(inplace = True, drop = True)
#	print(features.head())
	# drop all columns from domain file except the gene name and its domain
	domains = domains[["Gene", "Domain"]]
	# drop all rows that don't have a TE domain
	TE_domains = ["PF03184", "PF02914", "PF13358", "PF03732", "PF00665", "PF00077"]
	domains = domains[domains.Domain.isin(TE_domains)]
	# check that there are any genes with TE domains (if not, write the input gff to file without changing it)
	if domains.shape[0] == 0:
		features.to_csv("relabelled.gff", sep = "\t", index = False, header = False)
		print("No genes to relabel so original GFF file written to relabelled.gff")
		return 
	print(domains.shape[0], "genes with TE domains found")
#	print(domains.head())
	# collapse duplicate entries (which may be caused by split domain matches, so the same domain appears as >1 hit for the same gene)
	domains.drop_duplicates(inplace = True)
	# adjust gene names to strip off the _[0-9] added by interpro (NB right-hand split here should preserve underscores in gene names)
	new_gene_names = []
	for hit in domains.itertuples():
		new_gene_names.append(hit.Gene.rsplit("_", 1)[0])
	domains["Gene"] = new_gene_names
#	print(domains.head())
	# separate out entries for genes that have TE domains
	TEs_mislabelled_as_genes = features[features.Gene.isin(domains.Gene)]
	mislabelled_feature_count = TEs_mislabelled_as_genes.shape[0]
	# drop these entries from the original GFF
	features.drop(TEs_mislabelled_as_genes.index, inplace = True)
#	print("After screening out features for genes with TE domains, GFF has", features.shape[0], "entries")
	# add a column containing the TE domain to the mislabelled feature df
	TEs_mislabelled_as_genes = pd.merge(TEs_mislabelled_as_genes, domains)
	# check if each mislabelled feature has had only one domain added
	if mislabelled_feature_count == TEs_mislabelled_as_genes.shape[0]:
		print("Each mislabelled feature has had only one domain name added")
	elif mislabelled_feature_count > TEs_mislabelled_as_genes.shape[0]:
		print("ERROR: some features have not had a domain name added")
	else:
		print("ERROR: some features have had multiple domain names added")	
	# drop every feature that isn't a gene from the mislabelled features
	TEs_mislabelled_as_genes = TEs_mislabelled_as_genes[TEs_mislabelled_as_genes["Feature"] == "gene"]
#	print(TEs_mislabelled_as_genes)
	# relabel feature as TE, and name in the form: Chromosome;Start;End;TE_name;TE_domain;TE_family;TE_class
	TEs_mislabelled_as_genes["Feature"] = "TE"
	TEs_mislabelled_as_genes["Gene"] = TEs_mislabelled_as_genes[["Chromosome","Start","End","Gene","Domain"]].apply(lambda x: ";".join(x), axis = 1)
	TEs_mislabelled_as_genes["Gene"] = TEs_mislabelled_as_genes["Gene"] + ";Unknown;Unknown"
	TEs_mislabelled_as_genes.drop(["Domain"], inplace = True, axis = 1)
	# add relabelled TEs to the rest of the features and sort
	features = pd.concat([features, TEs_mislabelled_as_genes], ignore_index = True)
	features.reset_index(drop = True, inplace = True)
	features.sort_values(by = ["Chromosome", "Start", "Feature"], ascending = [True, True, False])
	# output to file
	features.to_csv("relabelled.gff", sep = "\t", index = False, header = False)
	print("Relabelled GFF file written to relabelled.gff")
	return

interpro_parse_gene(interpro = interpro_file, GFF = GFF_file)
