#!/usr/bin/env python

"""GenBank data ingestion into MySQL database. Ingests records directly from
    GenBank via a text file containing a list of accessions."""

## Imports ##
import argparse
import gb_csv_module as gcm

## Arguments ##
parser = argparse.ArgumentParser(
    description="Adding GenBank data to the database.")
req_group = parser.add_argument_group('required arguments')
req_group.add_argument('--db_user', dest='db_user', help="Database username",
                       required=True, metavar='USERNAME')
req_group.add_argument('--db_pass', dest='db_pass', help="Database password",
                       required=True, metavar='PASSWORD')

req_group.add_argument('-a', help="""Path to single-column text 
                    file containing accession numbers to search on GenBank.""",
                    dest='accessions', required=True)
req_group.add_argument('-e', help="""Enter your email address in order 
                    to access NCBI.""", dest='users_email', required=True)
req_group.add_argument('-p', help="""The prefix for the database id 
                    names.""", dest='prefix', required=True)
req_group.add_argument('-n', help="""The start number for the database 
                    id names""", dest='number', required=True)
req_group.add_argument('-z', help="""By how many zeros the number 
                    should be 0-padded.""", dest='padding', required=True)
parser.add_argument('-k', help="""Field of the genbank-file where the 
                    name of the entry can be found.""", dest='key',
                    default='LOCUS', choices=['LOCUS', 'ACCESSION',
                                              'DEFINITION'])

args = parser.parse_args()

## Functions ##

# Check login details
gcm.check_login_details(args.db_user, args.db_pass)

# Convert text-file of accessions to a list.
acc_list = gcm.text_to_list(args.accessions)

# Check the provided accessions are formatted correctly.
new_acc_list = gcm.check_acc_format(acc_list)

# Check accessions in list don't already exist in 'name' column of metadata
# db table or the 'rejected' table, dropping any that do.
filtered_accs = gcm.check_accs_in_db(new_acc_list)

# Take list of of IDs/accessions and return dict of corresponding GenBank
# records.
records = gcm.return_gb_data(filtered_accs, args.users_email)

# Check records for duplicates
filt_recs_1, reject_r1 = gcm.check_recs_for_dups(records)

# Check seqs don't already exist in the database
filt_recs_2, reject_r2 = gcm.check_seqs_in_db(filt_recs_1)

# Load rejected accessions into rejected_accs table
reject = reject_r1 + reject_r2
if reject:
    gcm.load_reject_table(reject)

# Create dict with old input ids (keys) and new database ids (values).
dict_new_ids = gcm.new_ids(filt_recs_2, args.prefix, args.number, args.padding)

# Check new ids are not already present in the database
gcm.check_ids(list(dict_new_ids.values()), 'ingest')

# Extract metadata from gb_records_dict and write to DataFrame.
gb_met_df = gcm.extract_metadata(filt_recs_2)

# In dataframe insert column with new database ids
gb_df_new_ids = gcm.change_names_csv(gb_met_df, dict_new_ids)

# Replace old input ID with new database ID in genbank file
gcm.change_ids_genbank(filt_recs_2, dict_new_ids, args.key)

# Add column with version number
gb_df_new_ids['version'] = 0

# Reorder dataframe columns for the database.
gb_df_reformatted = gcm.reformat_df_cols(gb_df_new_ids)

# Change the features for CDS
gcm.alter_features(filt_recs_2)

# Load new ids into master table
gcm.load_ids_to_master(list(dict_new_ids.values()))

# Push the genbank data into the database
gcm.load_gb_dict_into_db(filt_recs_2)

# Push the metadata into the database
gcm.load_df_into_db(gb_df_reformatted)

# Update master table
db_ids = list(dict_new_ids.values())

gcm.update_master_table(db_ids, db_ids, 'ingest')