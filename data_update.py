#!/usr/bin/env python

"""Data update in MySQL database. Adds new rows for each record and
    increments version number."""

## Imports ##
import argparse
import pandas as pd
import gb_csv_module as gcm
from argparse_formatter import FlexiFormatter

## Arguments ##
parser = argparse.ArgumentParser(description='Modifying data in the database',
                                 epilog="""
                                 --NOTES--
                                 * This script adds modified data pertaining to \
                                 existing records to the database. 
                                 * Genetic data and metadata can be updated \
                                 independently (i.e. one may parse a GenBank \
                                 file without a CSV or vice versa).
                                 * Version incrementation is automated in \
                                 this script, so the user should not manually \
                                 update version numbers when updating records \
                                 otherwise version numbers will be skipped.
                                 * Whatever records are ingested through this \
                                 script will be loaded into the database \
                                 regardless of whether they are duplicates of \
                                 pre-existing records (the database will not \
                                 recognise them as duplicates as all fields \
                                 with ‘unique’ constraints on them will have \
                                 been automatically incremented), so for \
                                 efficiency the user should ensure that any \
                                 unedited records are removed from the input \
                                 files before parsing them.""",
                                 formatter_class=FlexiFormatter)
req_group = parser.add_argument_group('required arguments')
req_group.add_argument('--db_user', dest='db_user', help="Database username",
                       metavar='', required=True)
req_group.add_argument('--db_pass', dest='db_pass', help="Database password",
                       metavar='', required=True)
parser.add_argument('-gb',
                    help="Path to genbank-file to ingest into the database.",
                    dest='input_genbank', metavar='')
parser.add_argument('-csv',
                    help="Path to genbank-file to ingest into the databse.",
                    dest='input_csv', metavar='')
parser.add_argument('-k',
                    help="""Field of the genbank-file where the name of the \
                    entry can be found.""", dest='key', default='LOCUS',
                    metavar='', choices=['LOCUS', 'ACCESSION', 'DEFINITION'])


#Subparser for optional manual update
#subparsers = parser.add_subparsers(dest="manual_update",
#                                   description="""Manually update records in
#                                   database.""")
#parser_manup = subparsers.add_parser('MANUP', help="""Manually update records
#                                    in database —— For help file see
#                                    'data_update.py MANUP -h'.""")
#parser_manup.add_argument('-s', '--specs', help="""Specifications of records
#                            to be updated. (E.g. 'length>25000'
#                            'country=United Kingdom') WARNING: IF THIS IS NOT
#                            PROVIDED, ALL DATABASE RECORDS WILL BE UPDATED.""",
#                          dest='mysql_specs', metavar='{specifications}',
#                          nargs='+')
#parser_manup.add_argument('-u', '--update', help="""Specifications of updates
#                            to be made. (E.g. 'locomotion=arboreal'
#                            'size=12mm')""", dest='update_specs',
#                          metavar='{update_specs}', nargs='+')
#parser_manup.add_argument('-q', '--custom_query', help="""Custom MySQL query to
#                            update data in database. (e.g. \"UPDATE metadata
#                            SET subregion='Sabah', size='12mm' WHERE
#                            name='BIOD00234';\") NOTE: Your custom specification
#                             must be enclosed by double quotations, as
#                             above.""", dest='custom_query',
#                          metavar='{custom_query}')

# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', "-gb", "./testdata/replace.gb", "-csv", "./testdata/replace.csv", '-k', 'LOCUS'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', "-gb", "./testdata/CHINAS.gb", "-csv", "./testdata/CHINAS.csv", '-k', 'LOCUS'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', "-csv", "./testdata/Honduras_testrun.csv", '-gb', './testdata/Honduras_testrun.gb', '-k', 'LOCUS'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', "MANUP", '-s', 'db_id=TEST099', '-u', 'size=;)'])

args = parser.parse_args()

"""
#This is for manually updating specific cells in database tables. 

if args.manual_update:
    # Define restrictions
    if args.custom_query and (args.mysql_specs or args.update_specs):
        parser.error("the following incompatible options are selected: -q, "
                     "[-s/-u]")
    if not args.custom_query and not args.update_specs:
        parser.error("the following argument is required: -u")

    if args.custom_query:

        mysql_query = args.custom_query

    else:

        mysql_query = gcm.construct_sql_update_query(None, args.update_specs, 
                                                     args.mysql_specs)
    #Case where the update and the spec share the same column?

    print(mysql_query)

    #Connect to database and execute command
    gcm.execute_query(mysql_query, args.db_user, args.db_pass)

    gcm.update_master_table(None, None, 'replace')
"""

## Functions ##

# Define restictions
if not (args.input_genbank or args.input_csv):
    parser.error("no update files provided.")

# Check login details
gcm.check_login_details(args.db_user, args.db_pass)

if args.input_genbank and args.input_csv:

    # Create dict from genbank-file
    gb_dict = gcm.gb_into_dictionary(args.input_genbank, args.key)

    # Create dataframe from csv file
    csv_df = pd.read_csv(args.input_csv, quotechar='"')

    # Check if the header in the csv is correct
    gcm.correct_header(csv_df, 'replace')

    # Check if latitude and longitude are formatted correctly
    gcm.lat_long(csv_df)

    # Check if the genbank and metadata file have matching entries
    csv_df, gb_dict = gcm.matching_inputids(csv_df, gb_dict, 'replace')

    # Generate ids list
    ids_list = list(gb_dict.keys())

    gb_ids, meta_ids = [ids_list, ids_list]

else:

    if args.input_genbank:

        csv_df = None

        # Create dict from genbank-file
        gb_dict = gcm.gb_into_dictionary(args.input_genbank, args.key)

        # Generate ids list
        ids_list = list(gb_dict.keys())

        gb_ids, meta_ids = [ids_list, None]

    if args.input_csv:

        gb_dict = None

        # Create dataframe from csv file
        csv_df = pd.read_csv(args.input_csv, quotechar='"')

        # Check if the header in the csv is correct
        gcm.correct_header(csv_df, 'replace')

        # Check if latitude and longitude are floats
        gcm.lat_long(csv_df)

        # Add version column
        csv_df['version'] = 0

        # Generate ids list
        ids_list = list(csv_df['db_id'])

        gb_ids, meta_ids = [None, ids_list]

# Check ids are already present in the database
gcm.check_ids(ids_list, 'replace')

# Update Data
gcm.update_data(csv_df, gb_dict)

# Update Master Table
gcm.update_master_table(gb_ids, meta_ids, 'update')

print('Done.')