#!/usr/bin/env python

"""Data ingestion into MySQL database.


Modifying and uploading a genbank-file and metadata in a csv-file into a MySQL database.
"""

import argparse
import pandas as pd
import sys

# Import the module:
import gb_csv_module as gcm


# Define arguments to be parsed
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

# args = parser.parse_args(["-gb", "/Users/lukeswaby-petts/Desktop/Work/Wildlife Research /Alfried/Mission 2/mtgendb/testdata/test_genbank.gb", "-csv", "/Users/lukeswaby-petts/Desktop/Work/Wildlife Research /Alfried/Mission 2/mtgendb/testdata/test_metadata.csv", "-r", "False", "-s", "Coleoptera", "-e", "luke.swaby@nhm.ac.uk", "-p", "TEST", "-n", "1", "-z", "3", "-k", "LOCUS"])


# Create dict from genbank-file and dataframe from csv-file

gb_dict = gcm.gb_into_dictionary(args.input_genbank, args.key)  #gb_dict = gcm.gb_into_dictionary(args.input_genbank, key="LOCUS")
#print("L: gb_into_dictionary() done")

csv_df = pd.read_csv(args.input_csv, quotechar = '"')               # PD.READ_CSV: Reads a comma-separated values (csv) file into DataFrame.
#print("L: read_csv() done")                                                                  # QUOTECHAR = "" - The character used to denote the start and end of a quoted item. Quoted items can include the delimiter and it will be ignored.

# Execute the module functions

#Check if the header in the csv is correct
gcm.correct_header(csv_df)
#print("L: correct_header() done")

#Check if the genbank and metadata file have matching entries
new_csv_df, new_gb_dict = gcm.matching_inputids(csv_df, gb_dict)
#print("L: matching_inputids() done (Well done Luke, ur a g)")
"""
Outputs new_csv_df and new_gb_dict
Need to somehow find a way to make subsequent funtions take these new files as their inputs.
"""

#Change the names of the IDs

#Create dictionary with old input ids (key) and new database ids (value)
dict_new_ids = gcm.new_ids(new_gb_dict, args.prefix, args.number, args.padding)    #dict_new_ids = gcm.new_ids(new_gb_dict, 'TEST', 1, 3)
#print("L: new_ids() done")

#Check if new ids are not already present in the database
# gcm.check_new_ids(dict_new_ids)

#In dataframe insert column with new database ids
df_new_ids = gcm.change_names_csv(new_csv_df, dict_new_ids)
#print("L: change_names_csv() done")

##Search for ncbi lineage with tax id, save custom lineages if not found on ncbi
lineages = gcm.get_ncbi_lineage(df_new_ids, args.users_email, args.searchterm)   # lineages = gcm.get_ncbi_lineage(df_new_ids, 'luke.swaby@nhm.ac.uk', 'Coleoptera')
#print("L: ncbi_lineage() done")

##User decides if he wants to reject entries with custom lineage information or not (-r True/False flag)

dict_accepted = gcm.rejecting_entries_new_gb(lineages, new_gb_dict, df_new_ids, args.reject_custom_lineage)
#print("L: rejecting_entries_new_gb() done")
df_accepted = gcm.rejecting_entries_new_csv(lineages, df_new_ids, args.reject_custom_lineage)
#print("L: rejecting_entries_new_csv() done")

#Create a new dictionary with the added taxonomy information
new_dict = gcm.insert_taxid(lineages, dict_accepted)
#print("L: insert_taxid() done")
#Replace old input ID with new database ID in genbank file
gcm.change_ids_genbank(new_dict, dict_new_ids, args.key)
#print("L: change_ids_genbank() done")
#Change the features for CDS
gcm.alter_features(new_dict)
#print("L: alter_features() done")

#Add columns with lineages and tax id to dataframe
df_with_lineages = gcm.add_lineage_df(df_accepted, lineages)
#print("L: df_with_lineages set to add_lineage_df(df_accepted, lineages)")

##Push the genbank data into the database

gcm.load_gb_dict_into_db(new_dict)
print("L: load_gb_dict_into_db() done")

##Push the metadata into the database

gcm.load_df_into_db(df_with_lineages)
print("L: load_df_into_db() done")


