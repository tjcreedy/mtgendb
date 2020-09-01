#!/usr/bin/env python

"""Data removal in MySQL database. Option to remove all records
    corresponding to input database IDs or just specified versions."""

## Imports ##
import argparse
import gb_csv_module as gcm

## Arguments ##
parser = argparse.ArgumentParser(description='Removing records in the database')
req_group = parser.add_argument_group('required arguments')
req_group.add_argument('--db_user', help="Database username", dest='db_user',
                       metavar='{db_username}', required=True)
req_group.add_argument('--db_pass', help="Database password", dest='db_pass',
                       metavar='{db_password}', required=True)

subparsers = parser.add_subparsers(dest='del_versions', description="""Choose 
                                    whether to delete all versions of a record 
                                    or just one.""")

# Parser for removal of all versions of target record(s)
remove_all = subparsers.add_parser('ALL', help="""Remove all versions of target 
                                    record(s) from all tables in database.""")
remove_all.add_argument('-t', '--txtfile', help="""Path to 1-column text file 
                        containing IDS of records to be deleted from 
                        database.""", dest='txtfile', metavar='{txtfile_path}')
remove_all.add_argument('-id', help="""Database ID(s) of record(s) to be 
                        deleted.""", dest='rec_id', metavar='{db_id}',
                        nargs='+')

# Parser for removal of single version of target record(s)
remove_version = subparsers.add_parser('VER', help="""Remove single version of 
                                        target record(s).""")
remove_version.add_argument('-t', '--txtfile', help="""Path to 3-column text 
                            file containing IDS and the versions of the records 
                            to be deleted from database - db_id <TAB> bio_ver 
                            <TAB> meta_ver.""", dest='txtfile',
                            metavar='{txtfile_path}')
remove_version.add_argument('-id', help="""Database ID of record to be 
                            deleted.""", dest='rec_id', metavar='{db_id}')
remove_version.add_argument('-m', '--met_ver', help="""Metadata version you wish 
                            to rollback to""", dest='meta_version',
                            metavar='{meta_version}', type=int)
remove_version.add_argument('-b', '--bio_ver', help="""Bioentry version to be 
                            rolled back to""", dest='bio_version',
                            metavar='{bio_version}', type=int)

#args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'VER', '-t', 'testdata/removes.txt'])

args = parser.parse_args()

## Functions ##

#Check login details
gcm.check_login_details(args.db_user, args.db_pass)

if args.del_versions == 'ALL':
    #Define restrictions
    if not (args.txtfile or args.rec_id):
        parser.error("one of the following arguments must be provided: -t/-id")

    db_ids = []

    if args.txtfile:

        db_ids = gcm.text_to_list(args.txtfile)

    if args.rec_id:

        db_ids.append(args.rec_id)

    gcm.check_ids(list(set(db_ids)), 'remove')

    gcm.remove_recs(db_ids)

else:
    #Define restrictions
    if not (args.txtfile or args.rec_id or args.version_no):
        parser.error("one of the following arguments must be provided: -t or "
                     "-id/-b/-m")
    if args.rec_id and not (args.bio_version or args.meta_version):
        parser.error('at least one of the following arguments is required: '
                     '-b/-m')
    if (args.bio_version or args.meta_version) and not args.rec_id:
        parser.error('the following argument is required: -id')

    versions_dict = {}

    if args.txtfile:
        #Convert to dict
        versions_dict = gcm.versions_to_dict(args.txtfile)

    if args.rec_id:

        versions_dict[args.rec_id] = {'b': args.bio_version,
                                      'm': args.meta_version}

    gcm.check_ids(list(set(versions_dict.keys())), 'remove')

    # Delete target records from database
    gcm.remove_versions(versions_dict)

print('Done.')

