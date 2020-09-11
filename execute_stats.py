#!/usr/bin/env python

"""Summary and statistics for genbank-files."""

## Imports ##
import gb_csv_module as gcm
import statistics as st
import pandas as pd
import subprocess
import argparse

## Arguments ##
parser = argparse.ArgumentParser(description="Summaries and Statistics")
parser.add_argument('-gb', help="""Path to GenBank file to ingest into the 
                    database.""", dest='input_genbank', metavar='',
                    required=True)
parser.add_argument('-csv', help="""Path to CSV file to ingest into the 
                    databse.""", dest='input_csv', metavar='', required=True)

args = parser.parse_args()

##Input gebank file as dict

gb_dict = gcm.gb_into_dictionary(args.input_genbank, key="LOCUS")

##Input csv-file as dataframe

csv_df = pd.read_csv(args.input_csv, quotechar='"')

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

gcm.alter_features(gb_dict)
contig_lengths = st.contig_length(gb_dict)

#Calculate mean, std, variance, median, range of gene length
st.length_variation(contig_lengths)

##Analyze the gene order
df_gene_order = st.gene_order(gb_dict)
order = st.extract_order(gb_dict)
st.analyze_order(order)

##Analyze family distribution
st.family_distribution(csv_df)

##Analyze country distribution
st.country_distribution(csv_df)

##Output plots generated in R
subprocess.call("Rscript create_plots.R", shell=True)

##Print gene sequences to fasta-file
st.gene_sequence(gb_dict)