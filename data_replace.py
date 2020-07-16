#!/usr/bin/env python

"""Data replacement in MySQL database.
"""

#Import modules
import argparse
import pandas as pd
import gb_csv_module as gcm

#Arguments to be parsed
parser = argparse.ArgumentParser(description='Modifying data in the database')
req_group = parser.add_argument_group('Required arguments')
req_group.add_argument('--db_user', dest='db_user', metavar='{db_username}', required=True, help="Database username")
req_group.add_argument('--db_pass', dest='db_pass', metavar='{db_password}', required=True, help="Database password")
parser.add_argument('-gb', '--genbankfile', dest='input_genbank', help="Name of genbank-file to ingest into the database.")
parser.add_argument('-csv', '--csvfile', dest='input_csv', help="Name of genbank-file to ingest into the databse.")
parser.add_argument('-k', '--key', dest='key', default='LOCUS', choices=['LOCUS', 'ACCESSION', 'DEFINITION'], help="Field of the genbank-file where the name of the entry can be found.")

#Subparser for optional manual update
subparsers = parser.add_subparsers(dest="manual_update", description='Manually update records in database.')
parser_manup = subparsers.add_parser('MANUP', help="Manually update records in database —— For help file see 'data_update.py MANUP -h'.")
parser_manup.add_argument('-s', '--specs', dest='mysql_specs', metavar='{specifications}', nargs='+', help="Specifications of records to be updated. (E.g. 'length>25000' 'country=United Kingdom') WARNING: IF THIS IS NOT PROVIDED, ALL DATABASE RECORDS WILL BE UPDATED.")
parser_manup.add_argument('-u', '--update', dest='update_specs', metavar='{update_specs}', nargs='+', help="Specifications of updates to be made. (E.g. 'locomotion=arboreal' 'size=12mm')")
parser_manup.add_argument('-q', '--custom_query', dest='custom_query', metavar='{custom_query}', help="Custom MySQL query to update data in database. (e.g. \"UPDATE metadata SET subregion='Sabah', size='12mm' WHERE name='BIOD00234';\") NOTE: Your custom specification must be enclosed by double quotations, as above.")

args = parser.parse_args()


#Define restrictions
if not (args.manual_update or args.input_genbank or args.input_csv):
    parser.error("no update option selected.")
if args.manual_update and (args.input_genbank or args.input_csv):
    parser.error("the following incompatible options are selected: MANUP, [-gb/-csv]")
if args.manual_update:
    if args.custom_query and (args.mysql_specs or args.update_specs):
        parser.error("the following incompatible options are selected: -q, [-s/-u]")
    if not args.custom_query and not args.update_specs:
        parser.error("the following argument is required: -u")

#Script
if args.manual_update:

    if args.custom_query:

        mysql_query = args.custom_query

    else:

        mysql_query = gcm.construct_sql_update_query(None, args.update, args.mysql_specs)
    #Case where the update and the spec share the same column?

    #Connect to database and execute command
    names_dict = gcm.execute_query(mysql_query, args.db_user, args.db_pass)

    # Add 1 to version no.?
#-----------
else:

    if args.input_genbank:

        #Create dict from genbank-file
        gb_dict = gcm.gb_into_dictionary(args.input_genbank, args.key)

        #Generate ids list
        ids_list = list(gb_dict.keys())

    if args.input_csv:

        #Create dataframe from csv file
        csv_df = pd.read_csv(args.input_csv, quotechar='"')

        #Check if the header in the csv is correct
        gcm.correct_header(csv_df)

        #Generate ids list
        ids_list = list(csv_df['name'])

    if args.input_genbank and args.input_csv:

        #Check if the genbank and metadata file have matching entries
        csv_df, gb_dict = gcm.matching_inputids(csv_df, gb_dict)

        #Generate ids list
        ids_list = list(gb_dict.keys())

    #Check if new ids are not already present in the database
    gcm.check_ids(args.db_user, args.db_pass, ids_list, 'replace')

    #Add 1 to version_no. of each new data using ids_list

"""
    if args.input_genbank:
    
        #Overwrite genetic info with gb_dict
        
        IT MUST:
        1. Overwrite data but assume same db_id (metadata.db_id, bioentry.name, bioentry.accession)
    
    if args.input_csv:
    
        #Overwrite metadata with csv_df
"""

mysql_query = "UPDATE table SET col1=newvalue1, col2=newvalue2 WHERE spec"
