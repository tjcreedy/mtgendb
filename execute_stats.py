#!/usr/bin/env python

"""Summary and statistics for genbank-files.
"""

import gb_csv_module as gcm
import statistics as st
import pandas as pd
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description="Summaries and Statistics")
parser.add_argument('-gb', '--genbankfile', dest='input_genbank', required=True, help="Name of GenBank file to ingest into the database.")
parser.add_argument('-csv', '--csvfile', dest='input_csv', required=True, help="Name of CSV file to ingest into the databse.")
args = parser.parse_args()

# args = parser.parse_args(["-gb", "/Users/lukeswaby-petts/Desktop/Work/Wildlife Research /Alfried/Mission 2/Testing/test_genbank.gb", "-csv", "/Users/lukeswaby-petts/Desktop/Work/Wildlife Research /Alfried/Mission 2/Testing/test_metadata_copy.csv"])

##Input gebank file as dict

gb_dict = gcm.gb_into_dictionary(args.input_genbank, key="LOCUS")
print("L: gb_into_dictionary() done")


##Input csv-file as dataframe

csv_df = pd.read_csv(args.input_csv, quotechar = '"')
print("L: read_csv() done")

##Extract features, topologies and sequence lengths for each entry

features = st.gb_features(gb_dict)
print("L: gb_features() done")

topologies = st.gb_topology(gb_dict)
print("L: gb_topology() done")

lengths = st.gb_length(gb_dict)
print("L: gb_length() done")


#Summaries in a dataframe
df_features = st.dict_to_df(lengths, topologies, features)
print("L: dict_to_df() done")

##Calculations on sequence lengths, topologies and features
#Ouput summary percentages

st.calculations(df_features)
print("L: calculations() done")
#Output mean, median, variance, std, range
st.calculations_ftr(df_features)
print("L: calculations_ftr() done")

##For each entry extract sequence length of each gene

gcm.alter_features(gb_dict) #Remove ambiguous gene names
print("L: alter_features() done")

contig_lengths = st.contig_length(gb_dict)
print("L: contig_length() done")

#Calculate mean, std, variance, median, range of gene length
st.length_variation(contig_lengths)
print("L: length_variation() done")

##Analyze the gene order
df_gene_order = st.gene_order(gb_dict)
print("L: gene_order() done")
order = st.extract_order(gb_dict)
print("L: extract_order() done")
st.analyze_order(order)
print("L: analyze_order() done")


# TITS UP FROM HERE

##Analyze family distribution
st.family_distribution(csv_df)
print("L: family_distribution() done")

##Analyze country distribution
st.country_distribution(csv_df)
print("L: country_distribution() done")

##Output plots generated in R
subprocess.call("Rscript EDITED_create_plots.r", shell=True)
print("L:  done")

##Print gene sequences to fasta-file
st.gene_sequence(gb_dict)
print("L: gene_sequence() done")