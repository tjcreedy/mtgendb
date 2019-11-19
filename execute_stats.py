#!/usr/bin/env python

"""Summary and statistics for genbank-files.
"""

import gb_csv_module as gcm
import statistics as st
import pandas as pd
import subprocess
import argparse

parser = argparse.ArgumentParser(description = "Summaries and Statistics")
parser.add_argument('-gb', '--genbankfile', dest = 'input_genbank', required = True, help = "Name of genbank-file to ingest into the database.")
parser.add_argument('-csv', '--csvfile', dest = 'input_csv', required = True, help = "Name of genbank-file to ingest into the databse.")
args = parser.parse_args()

##Input gebank file as dict

key = "LOCUS"
gb_dict = gcm.gb_into_dictionary(args.input_genbank, key)


##Input csv-file as dataframe

csv_df = pd.read_csv(args.input_csv, quotechar = '"')

##Extract features, topologies and sequence lengths for each entry

features = st.gb_features(gb_dict)
topologies = st.gb_topology(gb_dict)
lengths = st.gb_length(gb_dict)


#Summaries in a dataframe
df_features = st.dict_to_df(lengths, topologies, features)


##Calculations on sequence lengths, topologies and features
#Ouput summary percentages

st.calculations(df_features)

#Output mean, median, variance, std, range
st.calculations_ftr(df_features)


##For each entry extract sequence length of each gene
gcm.alter_features(gb_dict) #Remove ambiguous gene names
gene_lengths = st.gene_length(gb_dict)

#Calculate mean, std, variance, median, range of gene length
st.length_variation(gene_lengths)


##Analyze the gene order
df_gene_order = st.gene_order(gb_dict)
order = st.extract_order(gb_dict)
st.analyze_order(order)


##Analyze family distribution
st.family_distribution(csv_df)


##Analyze country distribution
st.country_distribution(csv_df)


##Output plots generated in R
subprocess.call("/usr/bin/Rscript create_plots.r", shell = True)


##Print gene sequences to fasta-file
st.gene_sequence(gb_dict)
