#!/usr/bin/env python


import argparse
import gb_csv_module as gcm

parser = argparse.ArgumentParser(description="Removing or rolling back records in the database")
subparsers = parser.add_subparsers(dest='action', description='Choose either data rollback or removal')

rollback = subparsers.add_parser('ROLLBACK', help='dshgksdkjhgds')
rollback.add_argument('-s', '--specs', dest='mysql_specs', default=[], metavar='{Specification(s)}', nargs='+', help="""Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') REMEMBER: Each individual spec must be enclosed by quotations and separated from the next by a space as above.""")

remove = subparsers.add_parser('REMOVE', help='dskhjfd')
remove.add_argument('-s', '--specs', dest='mysql_specs', default=[], metavar='{Specification(s)}', nargs='+', help="""Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') REMEMBER: Each individual spec must be enclosed by quotations and separated from the next by a space as above.""")

args = parser.parse_args()

#Constuct command
sql = gcm.construct_sql_output_query(None, ['name', 'db_id'], args.mysql_specs)

#Fetch names
names_dict = gcm.fetch_names(sql, args.db_user, args.db_pass)

#Record current meta and bio versions for each
versions = {db_id: {'metadata_version': gcm.check_latest_version(db_id)[0],
                    'bioentry_version': gcm.check_latest_version(db_id)[1]
                    } for db_id in names_dict.values()}

#Record minimum meta and bio version
min_bio_version = min([versions[name]['bioentry_version'] for name in versions.keys()])
min_meta_version = min([versions[name]['metadata_version'] for name in versions.keys()])

#Rollback bioentry

#Rollback metaata