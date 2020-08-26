#!/usr/bin/env python

"""Data rollback in MySQL database. Updates master table to highlight specified
    versions of records for data output."""

import argparse
import gb_csv_module as gcm

parser = argparse.ArgumentParser(description="Rolling back records in the database")
req_group = parser.add_argument_group('required arguments')
req_group.add_argument('--db_user', dest='db_user', metavar='{db_username}', required=True, help="Database username")
req_group.add_argument('--db_pass', dest='db_pass', metavar='{db_password}', required=True, help="Database password")
subparsers = parser.add_subparsers(dest='action', description='Choose either data rollback or removal')

single = subparsers.add_parser('SINGLE', help='Rollback data for single record specified on command line')
single.add_argument('-id', '--db_id', dest='db_id', metavar='{db_id', help='Database ID of record you wish to rollback.')
single.add_argument('-m', '--met_ver', dest='meta_version', metavar='{meta_version}', type=int, help='Metadata version you wish to rollback to')
single.add_argument('-b', '--bio_ver', dest='bio_version', metavar='{bio_version}', type=int, help='Bioentry version to be rolled back to')

multiple = subparsers.add_parser('MULTIPLE', help='Rollback data for multiple records specified in 3-column text file.')
multiple.add_argument('-t', '--txt', dest='text_file', metavar='{txtfile_path}', help='Path to 3-column text file: db_id <TAB> bio_ver <TAB> meta_ver')

args = parser.parse_args()

#Create target versions dict
if args.action == 'SINGLE':

    versions_dict = {args.db_id: {'b': args.bio_version, 'm': args.meta_version}}

else:

    versions_dict = gcm.versions_to_dict(args.text_file)

#Check ids
gcm.check_ids(args.db_user, args.db_pass, list(versions_dict.keys()), 'rollback')

#Update master table
gcm.rollback_versions(versions_dict)

print('Done.')