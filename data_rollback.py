#!/usr/bin/env python

"""Data rollback in MySQL database. Updates master table to highlight specified
    versions of records for data output."""

## Imports ##
import argparse
import gb_csv_module as gcm
from argparse_formatter import FlexiFormatter

## Arguments ##
parser = argparse.ArgumentParser(description="""Rolling back record versions \
                                in the database.""",
                                 epilog="""
                                 Each subparser contains its own arguments \
                                 and help file.\n""")
req_group = parser.add_argument_group('required arguments')
req_group.add_argument('--db_user', dest='db_user', help="Database username",
                       metavar='', required=True)
req_group.add_argument('--db_pass', dest='db_pass', help="Database password",
                       metavar='', required=True)

subparsers = parser.add_subparsers(dest='action', description="""Data rollback\
                                    for a single record on the command line or\
                                    multiple records in a text-file.""")

# Parser for rollback of single record on command line
single = subparsers.add_parser('SINGLE', help="""Rollback data for single \
                                record specified on command line.""")
single.add_argument('-id', help="""Database ID of record you wish to \
                    rollback.""", dest='db_id', metavar='')
single.add_argument('-m', help="""Metadata version number to be rolled back \
                    to""", dest='meta_version', metavar='', type=int)
single.add_argument('-b', help="""Bioentry version number to be rolled \
                    back to""", dest='bio_version', metavar='', type=int)

# Parser for rollback of multiple records by text file
multiple = subparsers.add_parser('MULTI', help="""Rollback data for multiple \
                                records specified in 3-column text file.""")
multiple.add_argument('-t', help="""Path to 3-column text file. Each row should\
                        consist of a tab-delimited list of the database ID of \
                        the target record, the bioentry version number to be \
                        rolled back to, and the metadata version number to be \
                        rolled back to:  db_id <TAB> bio_ver <TAB> meta_ver""",
                      dest='text_file',
                      metavar='')

args = parser.parse_args()

## Functions ##

#Check login details
gcm.check_login_details(args.db_user, args.db_pass)

#Create target versions dict
if args.action == 'SINGLE':
    versions_dict = {args.db_id: {'b': args.bio_version,
                                  'm': args.meta_version}}
else:
    versions_dict = gcm.versions_to_dict(args.text_file)

#Check ids
gcm.check_ids(list(versions_dict.keys()), 'rollback')

#Update master table
gcm.rollback_versions(versions_dict)

print('Done.')