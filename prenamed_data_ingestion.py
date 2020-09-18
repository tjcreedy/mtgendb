#!/usr/bin/env python

"""Data ingestion into MySQL database. Accepts a GenBank file containing
    genetic data and a csv containing metadata. Unlike the 'data_ingestion.py'
    script, this retains the original names of the input records."""

## Imports ##
import argparse
import pandas as pd
import gb_csv_module as gcm

## Arguments ##
parser = argparse.ArgumentParser(
    description="Adding a GenBank and a metadata CSV file to the database.")
req_group = parser.add_argument_group('required arguments')

req_group.add_argument('--db_user', dest='db_user', help="Database username",
                       required=True, metavar='USERNAME')
req_group.add_argument('--db_pass', dest='db_pass', help="Database password",
                       required=True, metavar='PASSWORD')
req_group.add_argument('-gb', help="""Name of genbank-file to ingest into 
                        the database.""", dest='input_genbank', required=True)
req_group.add_argument('-csv', help="""Name of genbank-file to ingest 
                    into the databse.""", dest='input_csv', required=True)
parser.add_argument('-r', help="""Set to 'False' if entries with 
                    custom lineage information should not be rejected.""",
                    dest='reject_custom_lineage', default='True',
                    choices=['True', 'False'], metavar='REJECT [True/False]')
parser.add_argument('-s', help="""Provide a taxonomic searchterm
                    to search for tax ids on NCBI for any entries without 
                    taxonomic data.""", dest='searchterm')
req_group.add_argument('-c', help="Path to ncbitaxid cache json file",
                       required=True, dest='taxidcache')
req_group.add_argument('-e', help="""Enter your email address in order 
                    to access NCBI.""", dest='users_email', required=True)
parser.add_argument('-k', help="""Field of the genbank-file where the 
                    name of the entry can be found.""", dest='key',
                    default='LOCUS', choices=['LOCUS', 'ACCESSION',
                                              'DEFINITION'])
args = parser.parse_args()

#args = parser.parse_args('--db_user root --db_pass mmgdatabase -gb testdata/gbmaster_2020-09-11_current.gb -csv testdata/AllMitogenomesMaster_noCERN_2020-09-17.csv -s Coleoptera -e thomas.creedy@gmail.com -c testdata/ncbitaxidcache.json'.split(' '))

## Functions ##

# Check login details
gcm.check_login_details(args.db_user, args.db_pass)

# Create dict from genbank-file and dataframe from csv-file
gb_dict = gcm.gb_into_dictionary(args.input_genbank, args.key)

csv_df = pd.read_csv(args.input_csv, quotechar='"')

# Check if the header in the csv is correct
gcm.correct_header(csv_df, 'ghost_ingest')

# Check if the genbank and metadata file have matching entries
new_csv_df, new_gb_dict = gcm.matching_inputids(csv_df, gb_dict, 'ghost_ingest')

# Check IDs don't already exist in database
db_ids = {rec.name: rec.name for rec in new_gb_dict.values()}
gcm.check_ids(list(db_ids.values()), 'ingest')

##Search for ncbi lineage with tax id, save custom lineages if not found on ncbi
lineages = gcm.get_ncbi_lineage(new_csv_df, args.taxidcache, 
                                args.users_email, args.searchterm)

#User decides whether to reject entries with custom lineage information or not
dict_accepted, df_accepted = gcm.rejecting_entries(
    lineages, new_gb_dict, new_csv_df, args.reject_custom_lineage,
    'ghost_ingest')

#Create a new dictionary with the added taxonomy information
new_dict = gcm.insert_taxid(lineages, dict_accepted)

#Replace old input ID with new database ID in genbank file
db_ids = {rec.name: rec.name for rec in new_dict.values()}
gcm.change_ids_genbank(new_dict, db_ids, args.key)

#Change the features for CDS
gcm.alter_features(new_dict)

#Add columns with lineages and tax id to dataframe
df_with_lineages = gcm.add_lineage_df(df_accepted, lineages)

#Add column with version number
df_with_lineages['version'] = 0

#Reorder metadata df columns to load into db.
df_with_lineages = gcm.reformat_df_cols(df_with_lineages)

#Load new ids into master table
db_ids = list(df_with_lineages['db_id']) #TODO: DELETE THIS LINE AND RETEST
gcm.load_ids_to_master(db_ids)

##Push the metadata into the database
gcm.load_df_into_db(df_with_lineages)

##Push the genbank data into the database
gcm.load_gb_dict_into_db(new_dict)

#Update master table
gcm.update_master_table(db_ids, db_ids, 'ingest')