#!/usr/bin/env python

"""GenBank data ingestion into MySQL database. Ingests records directly from
    GenBank via a text file containing a list of accessions."""

## Imports ##
import argparse
import gb_csv_module as gcm

## Arguments ##
parser = argparse.ArgumentParser(
    description="Adding GenBank data to the database.")

parser.add_argument('--db_user', help="Database username", dest='db_user',
                    metavar='db_username', required=True)
parser.add_argument('--db_pass', help="Database password", dest='db_pass',
                    metavar='db_password', required=True)
parser.add_argument('-a', '--accessions', help="""Path to text-file containing 
                    accession numbers to search on GenBank.""",
                    dest='input_accessions', required=True)
parser.add_argument('-e', '--email', help="""Enter your email address in order 
                    to access NCBI.""", dest='users_email', required=True)
parser.add_argument('-p', '--prefix', help="""The prefix for the database id 
                    names.""", dest='prefix', required=True)
parser.add_argument('-n', '--number', help="""The start number for the database 
                    id names""", dest='number', required=True)
parser.add_argument('-z', '--zeros', help="""By how many zeros the number 
                    should be 0-padded.""", dest='padding', required=True)
parser.add_argument('-k', '--key', help="""Field of the genbank-file where the 
                    name of the entry can be found.""", dest='key',
                    default='LOCUS', choices=['LOCUS', 'ACCESSION',
                                              'DEFINITION'])

args = parser.parse_args()

#LOCAL: args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', "-a", "/Users/lukeswaby-petts/Desktop/Work/Wildlife Research /Alfried/Mission 2/mtgendb/testdata/subset.txt", "-e", "luke.swaby@nhm.ac.uk", "-p", "GB", "-n", "1", "-z", "3", "-k", "LOCUS"])
# Long: args = parser.parse_args(['--db_user', 'root', 'db_pass', 'mmgdatabase', "-a", "/Users/lukeswaby-petts/Desktop/Work/Wildlife Research /Alfried/Mission 2/mtgendb/testdata/500accessions.txt", "-e", "luke.swaby@nhm.ac.uk", "-p", "GB", "-n", "1", "-z", "4", "-k", "LOCUS"])

#SERVER: args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', "-a", "/home/luke/Testing/subset.txt", "-e", "luke.swaby@nhm.ac.uk", "-p", "GB", "-n", "1", "-z", "3", "-k", "LOCUS"])
#COM: python3 gb_data_ingestion.py -a testdata/500accessions.txt -e luke.swaby@ngm.ac.uk -p GB -n 1 -z 3 -k LOCUS

## Functions ##

# Check login details
gcm.check_login_details(args.db_user, args.db_pass)

#Convert text-file of accessions to a list.
acc_list = gcm.text_to_list(args.input_accessions)

#Check the provided accessions are formatted correctly.
new_acc_list = gcm.check_acc_format(acc_list)

#Check accessions in list don't already exist in 'name' column of metadata
# db table, dropping any that do.
filtered_accs = gcm.check_accs_in_db(new_acc_list)

#Take list of of IDs/accessions and return dict of corresponding GenBank
# records.
records = gcm.return_gb_data(filtered_accs, args.users_email)

#Check records for duplicates
filt_recs_1, reject_r1 = gcm.check_recs_for_dups(records)

# Check seqs don't already exist in the database
filt_recs_2, reject_r2 = gcm.check_seqs_in_db(filt_recs_1)

#Load rejected accessions into rejected_accs table
reject = reject_r1 + reject_r2
if reject:
    gcm.load_reject_table(reject)

#Create dict with old input ids (keys) and new database ids (values).
dict_new_ids = gcm.new_ids(filt_recs_2, args.prefix, args.number, args.padding)

#Check new ids are not already present in the database
gcm.check_ids(list(dict_new_ids.values()), 'ingest')

#Extract metadata from gb_records_dict and write to DataFrame.
gb_met_df = gcm.extract_metadata(filt_recs_2)

#In dataframe insert column with new database ids
gb_df_new_ids = gcm.change_names_csv(gb_met_df, dict_new_ids)

#Replace old input ID with new database ID in genbank file
gcm.change_ids_genbank(filt_recs_2, dict_new_ids, args.key)

#Add column with version number
gb_df_new_ids['version'] = 0

#Reorder dataframe columns for the database.
gb_df_reformatted = gcm.reformat_df_cols(gb_df_new_ids)

#Change the features for CDS
gcm.alter_features(filt_recs_2)

#Load new ids into master table
gcm.load_ids_to_master(list(dict_new_ids.values()))

#Push the genbank data into the database
gcm.load_gb_dict_into_db(filt_recs_2)

#Push the metadata into the database
gcm.load_df_into_db(gb_df_reformatted)

#Update master table
db_ids = list(dict_new_ids.values())

gcm.update_master_table(db_ids, db_ids, 'ingest')