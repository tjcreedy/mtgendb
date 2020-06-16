#!/usr/bin/env python

"""Data retreival from MySQL database.
"""
import argparse
import pandas as pd
import sys

# Import the module:
import gb_csv_module as gcm


# Define arguments to be parsed
parser = argparse.ArgumentParser(description="Adding a genbank- and a csv-file to the database.")

#parser.add_argument('-sqlu', '--sqluser', dest = 'sql_user', required = True, help = "SQL username.")
#parser.add_argument('-sqlpw', '--sqlpassword', dest = 'sql_password', required = True, help = "SQL password.")
#parser.add_argument('-db', '--dbname', dest = 'db_name', required = True, help = "Database name.")

#parser.add_argument('-f', '--format', dest = 'output_format', required=True, choices=['genbank', 'fasta', 'csv'], help = "Output format")
#parser.add_argument('-t', '--table', dest = 'database_table', help = "Name of database table you wish to extract data from.")
#parser.add_argument('-c', '--columns', dest = 'table_columns', nargs='+', default = ['metadata.name', 'db_id'], help = "Name of table columns you wish to extract data from (e.g. name, .")
parser.add_argument('-o', '--out', dest = 'output_name', required = True, help = "Preferred filename for the output (extension will be added automatically according to your output format choice).")
#parser.add_argument('-s', '--specs', dest = 'mysql_specs', nargs='+', help = "Comma-separated list of mysql specifications (e.g. subregion='Sabah',collectionmethod='MALAISE').")
#parser.add_argument('-s', '--specs', dest = 'mysql_specs', type=str, help = "Comma-separated list of mysql specifications (e.g. subregion='Sabah',collectionmethod='MALAISE').")
subparsers = parser.add_subparsers(dest="subcommand")

parser_csv = subparsers.add_parser('CSV')
parser_csv.add_argument('-t', '--table', dest = 'database_table', choices=["metadata", "bioentry", "bioentry_dbxref", "bioentry_qualifier_value", "bioentry_reference", "biosequence", "comment", "taxon", "taxon_name"], help = "Name of database table you wish to extract data from.")
parser_csv.add_argument('-c', '--columns', dest = 'table_columns', nargs='+', default = '*', help = "Name of table columns you wish to extract data from (e.g. name, .")
parser_csv.add_argument('-csvs', '--csv_specs', dest = 'mysql_specs', nargs='+', help = "Comma-separated list of mysql specifications (e.g. subregion='Sabah',collectionmethod='MALAISE').")

parser_fasta = subparsers.add_parser('FA')
parser_fasta.add_argument('-fas', '--fasta_specifications', dest = 'mysql_specs', required=True, nargs='+', help = "Comma-separated list of mysql specifications (e.g. subregion='Sabah',collectionmethod='MALAISE').")

parser_gb = subparsers.add_parser('GB')
parser_gb.add_argument('-gbs', '--gb_specifications', dest = 'mysql_specs', required=True, nargs='+', help = "Comma-separated list of mysql specifications (e.g. subregion='Sabah',collectionmethod='MALAISE').")


"default = ['metadata.name', 'db_id']"

# args = parser.parse_args(["-sqlu", "root", "-sqlpw", "mmgdatabase", "-db", "mmg_test", "-t", "metadata", "-o", "metadateru", "-s", "subregion='Sabah'"])
# args = parser.parse_args(["-sqlu", "root", "-sqlpw", "mmgdatabase", "-db", "mmg_test", "-t", "metadata", "-c", "name", "length", "accession", "seq", "-o", "metadateru", "-s", "country='United Kingdon' description='Lucanus sp. BMNH 1425267 mitochondrion, complete genome'", "-f", "csv"])

args=parser.parse_args()

#gcm.csv_from_sql(args.sql_user, args.sql_password, args.db_name, args.database_table, args.table_columns, args.mysql_specs, args.output_name)
"""
print(f"OUT_NAME: {args.output_name}")
if args.subcommand == 'CSV':
    print(f"TABLE: {args.database_table}")
    print(f"COLUMNS: {args.table_columns}")
    print(f"SPECS: {args.mysql_specs}")
elif args.subcommand == 'FA':
    print(f"SPECS: {args.mysql_specs}")
elif args.subcommand == 'GB':
    print(f"SPECS: {args.mysql_specs}")
else:
    print("No output format specified.")
"""


if args.subcommand == 'CSV':

    mysql_command = gcm.construct_sql_command(args.database_table, args.table_columns, args.mysql_specs)

    gcm.csv_from_sql(mysql_command, args.output_name)

    # table, cols, spec = [None, '*', ['country=United Kingdom']]
    # table, cols, spec = ['metadata', '*', ['country=United Kingdom', 'description=Lucanus sp. BMNH 1425267 mitochondrion, complete genome']]
    # table, cols, spec = ['metadata', ['name', 'taxon_id'], ['country=United Kingdom']]
    # table, cols, spec = ['metadata', ['name', 'taxon_id'], None]
    # table, cols, spec = [None, ['name', 'taxon_id'], ['country=United Kingdom', 'description=Lucanus sp. BMNH 1425267 mitochondrion, complete genome']]
    # table, cols, spec = ['metadata', '*', ['country=United Kingdom', 'description=Lucanus sp. BMNH 1425267 mitochondrion, complete genome']]

    #table, cols, spec = ['metadata', ['name', 'latitude'], ['country=United Kingdom', 'description=Lucanus sp. BMNH 1425267 mitochondrion, complete genome']]
    # metadata.name, metadata.latitude ... ... WHERE (metadata.country= ...) AND (bioentry.description= ...);

else:

    mysql_command = gcm.construct_sql_command(None, ['name', 'db_id'], args.mysql_specs)


    # table, cols, spec = [None, ['metadata.name', 'db_id'], ['country=United Kingdom']]

    # table, cols, spec = [None, ['metadata.name', 'db_id'], ['country=United Kingdom', 'description=Lucanus sp. BMNH 1425267 mitochondrion, complete genome']]

    names_dict = gcm.fetch_names(mysql_command)

    records = gcm.fetch_recs(names_dict)

    gcm.seqfile_from_sql(records, args.output_name, args.output_format)

#print(f"FORMAT: {args.output_format}")

print(mysql_command)



"""
if ((args.database_table is not None) or (args.table_columns is not None)) and args.output_format != 'csv':

    x = input("ERROR: If table (-t) or columns (-c) arguments are provided, the script will assume you want a CSV. If you want anything else, then please provide just the specifications arg (-s).\nWould you like to switch output format to CSV? 'y'/'n'\n?>").lower()
    while not (x == 'y' or x == 'n'):
        x = input("ERROR: Please type 'y' to change output file format to CSV or 'n' to cancel the operation.").lower()
    if x == 'y':
        args.output_format = 'CSV'
    else:
        parser.error("Operation cancelled.")
"""

"""
INPUTS: CSV -t metadata                                 ---> ['metadata']                    ---> SELECT * FROM metadata;
        CSV -t metadata -c name db_id taxon_id       ---> ['metadata'] ["name", "db_id", "taxon_id"]       ---> SELECT name, db_id, taxon_id FROM metadata;
        CSV -c name length                                  ---> SELECT name, length FROM metadata JOIN bioentry .... =bioentry.bioentry_id;
        CSV -c name length -csvs country='United Kingdom'   ---> SELECT name, length FROM metadata JOIN bioentry .... =bioentry.bioentry_id WHERE metadata.coountry='United Kingdom';
        
        FA -fas country='United Kingdom'
        FA -fas country='United Kingdom' description='Lucanus sp. BMNH 1425267 mitochondrion, complete genome'
        
        GB -gbs country='United Kingdom'
        GB -gbs country='United Kingdom' description='Lucanus sp. BMNH 1425267 mitochondrion, complete genome'
         
"""


"""
mysql_command = gcm.construct_sql_command(args.table_columns, args.mysql_specs)


if args.subcommand == 'csv':

    gcm.csv_from_sql(mysql_command, args.output_name)

else:

    names_dict = gcm.fetch_names(mysql_command)

    records = gcm.fetch_recs(names_dict)

    gcm.seqfile_from_sql(records, args.output_name, args.output_format)
"""

"""
if FORMAT == CSV:
Need table and/or columns (and spec)

ELSE:
Need spec
    
    
"""




"""
TO-DO

1. What if specs are to do with numbers and not strings? (e.g. instead of country='UK', length>3000)
2. Which tables join which and how?
3. What happens if they have chosen multiple columns but have chosen an output format other than CSV?
4. If 'tables' provided, then set cols to '*' (as output will be CSV) [if table or cols given then format must be CSV]
5. 'name' column in both metadata and taxon_name. Must distinguish.
6. Sort out help page.

"""
