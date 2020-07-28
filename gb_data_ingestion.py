#!/usr/bin/env python

"""GenBank data ingestion into MySQL database.
"""

import argparse

# Import the module:
import gb_csv_module as gcm


# Define arguments to be parsed
parser = argparse.ArgumentParser(description="Adding GenBank data to the database.")

parser.add_argument('-a', '--accessions', dest='input_accessions', required=True, help="Name of text-file containing accession numbers to search on GenBank.")
parser.add_argument('-e', '--email', dest='users_email', required=True, help="Enter your email address in order to access NCBI.")
parser.add_argument('-p', '--prefix', dest='prefix', required=True, help="The prefix for the database id names.")
parser.add_argument('-n', '--number', dest='number', required=True, help="The start number for the database id names")
parser.add_argument('-z', '--zeros', dest='padding', required=True, help="By how many zeros the number should be 0-padded.")
parser.add_argument('-k', '--key', dest='key', default='LOCUS', choices=['LOCUS', 'ACCESSION', 'DEFINITION'], help="Field of the genbank-file where the name of the entry can be found.")

args = parser.parse_args()

#LOCAL: args = parser.parse_args(["-a", "/Users/lukeswaby-petts/Desktop/Work/Wildlife Research /Alfried/Mission 2/mtgendb/testdata/subset.txt", "-e", "luke.swaby@nhm.ac.uk", "-p", "GB", "-n", "1", "-z", "3", "-k", "LOCUS"])
# Long: args = parser.parse_args(["-a", "/Users/lukeswaby-petts/Desktop/Work/Wildlife Research /Alfried/Mission 2/mtgendb/testdata/500accessions.txt", "-e", "luke.swaby@nhm.ac.uk", "-p", "GB", "-n", "1", "-z", "4", "-k", "LOCUS"])

#SERVER: args = parser.parse_args(["-a", "/home/luke/Testing/subset.txt", "-e", "luke.swaby@nhm.ac.uk", "-p", "GB", "-n", "1", "-z", "3", "-k", "LOCUS"])
#COM: python3 gb_data_ingestion.py -a testdata/500accessions.txt -e luke.swaby@ngm.ac.uk -p GB -n 1 -z 3 -k LOCUS

#Convert text-file of accessions to a list.
acc_list = gcm.text_to_list(args.input_accessions)

#Check the provided accessions are formatted correctly.
new_acc_list = gcm.check_acc_format(acc_list)

#Take list of of IDs/accessions and return dict of corresponding GenBank records.
records = gcm.return_gb_data(new_acc_list, args.users_email)

#Create dict with old input ids (keys) and new database ids (values).
dict_new_ids = gcm.new_ids(records, args.prefix, args.number, args.padding)

#Extract metadata from gb_records_dict and write to DataFrame.
gb_met_df = gcm.extract_metadata(records)

#In dataframe insert column with new database ids
gb_df_new_ids = gcm.change_names_csv(gb_met_df, dict_new_ids)

#Replace old input ID with new database ID in genbank file
gcm.change_ids_genbank(records, dict_new_ids, args.key)

#Reorder dataframe columns for the database.
gb_df_reformatted = gcm.reformat_df_cols(gb_df_new_ids)

#Change the features for CDS
gcm.alter_features(records)

#Push the genbank data into the database
gcm.load_gb_dict_into_db(records)

#Push the metadata into the database
gcm.load_df_into_db(gb_df_reformatted)