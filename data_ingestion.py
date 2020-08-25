#!/usr/bin/env python

"""Data ingestion into MySQL database. Accepts a GenBank file containing
    genetic data and a csv containing metadata."""

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
parser.add_argument('-p', '--prefix', help="""The prefix for the database id 
                    names.""", dest='prefix', required=True)
parser.add_argument('-n', '--number', help="""The start number for the database 
                    id names""", dest='number', required=True)
parser.add_argument('-z', '--zeros', help="""By how many zeros the number should
                    be 0-padded.""", dest='padding', required=True)
parser.add_argument('-k', '--key', help="""Field of the genbank-file where the 
                    name of the entry can be found.""", dest='key',
                    default='LOCUS', choices=['LOCUS', 'ACCESSION',
                                              'DEFINITION'])

args = parser.parse_args()

# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', "-gb", "./testdata/test_genbank.gb", "-csv", "./testdata/test_metadata.csv", "-r", "False", "-s", "Coleoptera", "-e", "luke.swaby@nhm.ac.uk", "-p", "TEST", "-n", "1", "-z", "3", "-k", "LOCUS"])
#COM: python3 data_ingestion.py --db_user root --db_pass mmgdatabase -gb testdata/test_genbank.gb -csv testdata/test_metadata.csv -r False -s Coleoptera -e luke.swaby@nhm.ac.uk -p TEST -n 1 -z 3 -k LOCUS

## Functions ##

#Create dict from genbank-file and dataframe from csv-file
gb_dict = gcm.gb_into_dictionary(args.input_genbank, args.key)

csv_df = pd.read_csv(args.input_csv, quotechar='"')

#Check if the header in the csv is correct
gcm.correct_header(csv_df, 'ingest')

#Check if the genbank and metadata file have matching entries
new_csv_df, new_gb_dict = gcm.matching_inputids(csv_df, gb_dict, 'ingest')

#Create dictionary with old input ids (keys) and new database ids (values)
dict_new_ids = gcm.new_ids(new_gb_dict, args.prefix, args.number, args.padding)

#Check new ids are not already present in the database
gcm.check_ids(args.db_user, args.db_pass, list(dict_new_ids.values()), 'ingest')

#In dataframe insert column with new database ids
df_new_ids = gcm.change_names_csv(new_csv_df, dict_new_ids)

##Search for ncbi lineage with tax id, save custom lineages if not found on ncbi
lineages = gcm.get_ncbi_lineage(df_new_ids, args.users_email, args.searchterm)

#User decides whether to reject entries with custom lineage information or not
dict_accepted, df_accepted = gcm.rejecting_entries(
    lineages, new_gb_dict, df_new_ids, args.reject_custom_lineage)

#Create a new dictionary with the added taxonomy information
new_dict = gcm.insert_taxid(lineages, dict_accepted)

#Replace old input ID with new database ID in genbank file
gcm.change_ids_genbank(new_dict, dict_new_ids, args.key)

#Change the features for CDS
gcm.alter_features(new_dict)

#Add columns with lineages and tax id to dataframe
df_with_lineages = gcm.add_lineage_df(df_accepted, lineages)

#Add column with version number
df_with_lineages['version'] = 0

#Reorder metadata df columns to load into db.
df_with_lineages = gcm.reformat_df_cols(df_with_lineages)

#Load new ids into master table
gcm.load_ids_to_master(dict_new_ids)

##Push the metadata into the database
gcm.load_df_into_db(df_with_lineages)

##Push the genbank data into the database
gcm.load_gb_dict_into_db(new_dict)

#Update master table
db_ids = list(dict_new_ids.values())

gcm.update_master_table(db_ids, db_ids, 'ingest')