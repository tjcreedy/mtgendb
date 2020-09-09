#!/usr/bin/env python

"""Data ingestion into MySQL database. Accepts a GenBank file containing
    genetic data and a csv containing metadata. Unlike the 'data_ingestion.py'
    script, this retains the original names of the input records."""

#TODO: TEST THIS!

## Imports ##
import argparse
import pandas as pd
import gb_csv_module as gcm

## Arguments ##
parser = argparse.ArgumentParser(
    description="Adding a GenBank and a metadata CSV file to the database.")

parser.add_argument('--db_user', help="Database username", dest='db_user',
                    metavar='db_username', required=True)
parser.add_argument('--db_pass', help="Database password", dest='db_pass',
                    metavar='db_password', required=True)
parser.add_argument('-gb', '--genbankfile', help="""Name of genbank-file to 
                    ingest into the database.""", dest='input_genbank',
                    required=True)
parser.add_argument('-csv', '--csvfile', help="""Name of genbank-file to ingest 
                    into the databse.""", dest='input_csv', required=True)
parser.add_argument('-r', '--reject', help="""Set to 'false' if entries with 
                    custom lineage information should not be rejected.""",
                    dest='reject_custom_lineage', default='True',
                    choices=['True', 'False'])
parser.add_argument('-s', '--searchterm', help="""Provide a taxonomic searchterm
                    to search for tax ids on NCBI for any entries without 
                    taxonomic data (e.g. Coleoptera).""", dest='searchterm')
parser.add_argument('-e', '--email', help="""Enter your email address in order 
                    to access NCBI.""", dest='users_email', required=True)
parser.add_argument('-k', '--key', help="""Field of the genbank-file where the 
                    name of the entry can be found.""", dest='key',
                    default='LOCUS', choices=['LOCUS', 'ACCESSION',
                                              'DEFINITION'])

args = parser.parse_args()

## Functions ##

# Check login details
gcm.check_login_details(args.db_user, args.db_pass)

# Create dict from genbank-file and dataframe from csv-file
gb_dict = gcm.gb_into_dictionary(args.input_genbank, args.key)

csv_df = pd.read_csv(args.input_csv, quotechar='"')

# Check if the header in the csv is correct
gcm.correct_header(csv_df, 'ingest')

# Check if the genbank and metadata file have matching entries
new_csv_df, new_gb_dict = gcm.matching_inputids(csv_df, gb_dict, 'ingest')

# Check IDs don't already exist in database
db_ids = [rec.name for rec in gb_dict.values()]
gcm.check_ids(db_ids, 'ingest')

# In dataframe insert column with new database ids
new_csv_df['db_id'] = new_csv_df['name']

##Search for ncbi lineage with tax id, save custom lineages if not found on ncbi
lineages = gcm.get_ncbi_lineage(new_csv_df, args.users_email, args.searchterm)

#User decides whether to reject entries with custom lineage information or not
dict_accepted, df_accepted = gcm.rejecting_entries(
    lineages, new_gb_dict, new_csv_df, args.reject_custom_lineage)

#Create a new dictionary with the added taxonomy information
new_dict = gcm.insert_taxid(lineages, dict_accepted)

#Change the features for CDS
gcm.alter_features(new_dict)

#Add columns with lineages and tax id to dataframe
df_with_lineages = gcm.add_lineage_df(df_accepted, lineages)

#Add column with version number
df_with_lineages['version'] = 0

#Reorder metadata df columns to load into db.
df_with_lineages = gcm.reformat_df_cols(df_with_lineages)

#Load new ids into master table
gcm.load_ids_to_master(db_ids)

##Push the metadata into the database
gcm.load_df_into_db(df_with_lineages)

##Push the genbank data into the database
gcm.load_gb_dict_into_db(new_dict)

#Update master table
gcm.update_master_table(db_ids, db_ids, 'ingest')