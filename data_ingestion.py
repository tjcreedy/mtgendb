#!/usr/bin/env python

"""Data ingestion into MySQL database.

Modifying and uploading a genbank-file and metadata in a csv-file into a MySQL database.
"""

import argparse
import pandas as pd

##Import the module:
import gb_csv_module as gcm


##Define arguments to be parsed
parser = argparse.ArgumentParser(description="Adding a genbank- and a csv-file to the database.")

parser.add_argument('-gb', '--genbankfile', dest = 'input_genbank', required = True, help = "Name of genbank-file to ingest into the database.")
parser.add_argument('-csv', '--csvfile', dest = 'input_csv', required = True, help = "Name of genbank-file to ingest into the databse.")
parser.add_argument('-r', '--reject', dest = 'reject_custom_lineage', default = 'True', choices = ['True', 'False'], help = "Set to 'false' if entries with custom lineage information should not be rejected.")
parser.add_argument('-s', '--searchterm', dest = 'searchterm', help = "Provide a searchterm for the search for tax ids on NCBI.")
parser.add_argument('-e', '--email', dest = 'users_email', required = True, help = "Enter your email address in order to access NCBI.")
parser.add_argument('-p', '--prefix', dest = 'prefix', required = True, help = "The prefix for the database id names.")
parser.add_argument('-n', '--number', dest = 'number', required = True, help = "The start number for the database id names")
parser.add_argument('-z', '--zeros', dest = 'padding', required = True, help = "By how many zeros the number should be 0-padded.")
parser.add_argument('-k', '--key', dest = 'key', default = 'LOCUS', choices = ['LOCUS', 'ACCESSION', 'DEFINITION'], help = "Field of the genbank-file where the name of the entry can be found.")

args=parser.parse_args()


##Create dict from genbank-file and dataframe from csv-file

gb_dict = gcm.gb_into_dictionary(args.input_genbank, args.key)
csv_df = pd.read_csv(args.input_csv, quotechar = '"')


##Execute the module functions

#Check if the header in the csv is correct
gcm.correct_header(csv_df)

#Check if the genbank and metadata file have matching entries
gcm.matching_inputids(csv_df, gb_dict)

#Check if the metadata file has unique entry ids only
gcm.unique_csv_ids(csv_df)


##Change the names of the IDs

#Create dictionary with old input ids (key) and new database ids (value)
dict_new_ids = gcm.new_ids(gb_dict, args.prefix, args.number, args.padding)

#Check if new ids are not already present in the database
gcm.check_new_ids(dict_new_ids)

#In dataframe insert column with new database ids
df_new_ids = gcm.change_names_csv(csv_df, dict_new_ids)


##Search for ncbi lineage with tax id, save custom lineages if not found on ncbi

lineages = gcm.ncbi_lineage(df_new_ids, args.users_email, args.searchterm)


##User decides if he wants to reject entries with custom lineage information or not (-r True/False flag)

dict_accepted = gcm.rejecting_entries_new_gb(lineages, gb_dict, df_new_ids, args.reject_custom_lineage)
df_accepted = gcm.rejecting_entries_new_csv(lineages, df_new_ids, args.reject_custom_lineage)


#Create a new dictionary with the added taxonomy information
new_dict = gcm.insert_taxid(lineages, dict_accepted)

#Replace old input ID with new database ID in genbank file
gcm.change_ids_genbank(new_dict, dict_new_ids, args.key)

#Change the features for CDS
gcm.alter_features(new_dict)


#Add columns with lineages and tax id to dataframe
df_with_lineages = gcm.add_lineage_df(df_accepted, lineages)


##Push the genbank data into the database

gcm.load_gb_dict_into_db(new_dict)


##Push the metadata into the database

gcm.load_df_into_db(df_with_lineages)


