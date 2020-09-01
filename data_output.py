#!/usr/bin/env python

"""Data retreival from MySQL database. Outputs data in specified format
    according to SQL specifications."""

## Imports ##
import argparse
import re
import gb_csv_module as gcm


## Arguments ##
parser = argparse.ArgumentParser(
    description="Extracting data from the database.")
req_group = parser.add_argument_group('required arguments')
req_group.add_argument('--db_user', dest='db_user', help="Database username",
                       metavar='{db_username}', required=True)
req_group.add_argument('--db_pass', dest='db_pass', help="Database password",
                       metavar='{db_password}', required=True)
parser.add_argument('--all', help="""Use this flag if you wish to pull all 
                    versions of each record satisfying your query. By default 
                    this script will only pull current versions if flag is 
                    omitted.""", action='store_true')
parser.add_argument('-q', '--query', help="""Custom MySQL query to extract data 
                    from database. (e.g. \"SELECT * FROM metadata WHERE 
                    country='United Kingdom';\") NOTE: Your custom specification 
                    must be enclosed by double quotations, as above.""",
                    dest='mysql_query', metavar='custom_MySQL_query')

subparsers = parser.add_subparsers(dest="output_format",
                                   description='Desired output format:')

#Create the parser for the 'COUNT' command
parser_count = subparsers.add_parser('COUNT', help="""Prints an integer on the 
                                    command line —— For help file see 
                                    'data_output.py COUNT -h'.""")
parser_count.add_argument('-s', '--specs', help="""Comma-separated list of mysql
                            specifications.\n(e.g. 'subregion=Sabah' 
                            'length>25000' 'order=Coleoptera') REMEMBER: Each 
                            individual spec must be enclosed by quotations and 
                            separated from the next by a space as above.""",
                          dest='mysql_specs', default=[],
                          metavar='{Specification(s)}', nargs='+')
parser_count.add_argument('-x', '--taxonomy', help="""Taxonomic searchterm to 
                            run through database""", dest='taxonomy_spec',
                          metavar='{taxonomy}')

#Create the parser for the 'CSV' command
parser_csv = subparsers.add_parser('CSV', help="""Ouputs a .csv file —— For 
                                    help file see 'data_output.py CSV -h'.""")
parser_csv.add_argument('-o', '--out', help="""Preferred filename for the output
                        (extension will be added automatically according to 
                        your output format choice).""", dest='output_name',
                        metavar='{Output filename}', required=True)
parser_csv.add_argument('-t', '--table', help="""Name of database table you wish
                        to extract data from. (e.g. metadata.) If provided, 
                        the script will assume you require data from every 
                        column in the specified table. If this is not the case 
                        then stating only the required columns under the -c flag
                         will suffice.""", dest='database_table',
                        metavar='{Table name}',
                        choices=["metadata", "bioentry", "bioentry_dbxref",
                                 "bioentry_qualifier_value",
                                 "bioentry_reference", "biosequence",
                                 "seqfeature", "comment", "taxon", "taxon_name"])
parser_csv.add_argument('-c', '--columns', help="""Name of table columns you 
                        wish to extract data from. (e.g. name length 
                        description)""", dest='table_columns', default=['*'],
                        metavar='{Column name(s)}', nargs='+')
parser_csv.add_argument('-s', '--specs', help="""Comma-separated list of mysql 
                        specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 
                        'order=Coleoptera') NOTE: Each individual specification 
                        must be enclosed by quotations and separated from the 
                        next by a space as above.""", dest='mysql_specs',
                        default=[], metavar='{Specification(s)}', nargs='+')
parser_csv.add_argument('-x', '--taxonomy', help="""Taxonomic searchterm to run 
                        through database""", dest='taxonomy_spec',
                        metavar='{taxonomy}')

#Create the parser for the 'FASTA' command
parser_fasta = subparsers.add_parser('FASTA', help="""Ouputs a .fasta file —— 
                                    For help file see 'data_output.py FASTA 
                                    -h'.""")
parser_fasta.add_argument('-o', '--out', help="""Preferred filename for the
                            output (extension will be added automatically 
                            according to your output format choice).""",
                          dest='output_name', metavar='{Output filename}',
                          required=True)
parser_fasta.add_argument('-s', '--specs', help="""Comma-separated list of mysql
                            specifications.\n(e.g. 'subregion=Sabah' 
                            'length>25000' 'order=Coleoptera') NOTE: Each 
                            individual specification must be enclosed by 
                            quotations and separated from the next by a space as 
                            above.""", dest='mysql_specs', default=[],
                          metavar='{Specification(s)}', nargs='+')
parser_fasta.add_argument('-x', '--taxonomy', help="""Taxonomic searchterm to 
                            run through database""", dest='taxonomy_spec',
                          metavar='{taxonomy}')
parser_fasta.add_argument('-g', '--genes',  help="""Name of mitochondrial genes 
                            you wish to extract.', dest='genes""",
                          metavar='{gene_names}', nargs='+',
                          choices=['*', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3',
                                   'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L',
                                   'ND5', 'ND6'])

#Create the parser for the 'GB' command
parser_gb = subparsers.add_parser('GB', help="""Outputs an annotated .gb file 
                                —— For help file see 'data_output.py GB -h'.""")
parser_gb.add_argument('-o', '--out', help="""Preferred filename for the output 
                        (extension will be added automatically according to your 
                        output format choice).""", dest='output_name',
                       metavar='{Output filename}', required=True)
parser_gb.add_argument('-s', '--specs', help="""Comma-separated list of mysql 
                        specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 
                        'order=Coleoptera') NOTE: Each individual specification 
                        must be enclosed by quotations and separated from the 
                        next by a space as above.""", dest='mysql_specs',
                       default=[], metavar='{Specification(s)}', nargs='+')
parser_gb.add_argument('-x', '--taxonomy', help="""Taxonomic searchterm to run 
                        through database""", dest='taxonomy_spec',
                       metavar='{taxonomy}')
parser_gb.add_argument('-g', '--genes', help="""Name of mitochondrial genes you 
                        wish to extract.""", dest='genes',
                       metavar='{gene_names}', nargs='+',
                       choices=['*', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3',
                                'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L',
                                'ND5', 'ND6'])


# args = parser.parse_args(["--db_user", "root", "--db_pass", "mmgdatabase", '-q', "SELECT * FROM metadata WHERE country='China';", 'CSV', "-t", "metadata", '-c', 'name', "-o", "metadateru", "-s", "subregion='Sabah'"])
# args = parser.parse_args(["-sqlu", "root", "-sqlpw", "mmgdatabase", "-db", "mmg_test", "-t", "metadata", "-c", "name", "length", "accession", "seq", "-o", "metadateru", "-s", "country='United Kingdon' description='Lucanus sp. BMNH 1425267 mitochondrion, complete genome'", "-f", "csv"])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'CSV', '-o', 'Honduras_testrun', '-t', 'metadata', '-s', 'country=Honduras'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'CSV', '-o', 'gjhg', '-c', 'name', 'db_id', 'country', 'length', '-s', 'country=Honduras', 'length>4000'])

# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'CSV', '-o', 'Honduras_testrun', '-t', 'metadata', '-s', "db_id IN ('TEST001','TEST002')"])

# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'GB', '-o', 'finalesz007', '-x', 'Cucujiformia', '-s', 'country!=United Kingdom', '-g', 'COX1'])

# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'COUNT', '-x', 'Cucujiformia'])

# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'FASTA', '-o', 'outksis', '-s', 'country=Honduras', '-g', '*'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'FASTA', '-o', 'outksis', '-s',  'country!=Malaysia', '-g' , '*'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'FASTA', '-o', 'malays', '-s',  'species=Stenus boops', '-g', 'COX2', 'ND3', 'ATP6'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'CSV', '-o', 'Pan', '-t', 'metadata', '-s',  'country=Panama'])

args = parser.parse_args()

## Functions ##

#Define restrictions
if args.mysql_query and (args.taxonomy_spec or args.mysql_specs
                         or args.database_table or args.table_columns):
    parser.error("""query-building specifications and custom mysql query both 
    provided. Please proceed with one method.""")

#Check login details
gcm.check_login_details(args.db_user, args.db_pass)

if args.output_format == 'CSV':

    if args.mysql_query:

        mysql_command = args.mysql_query

    else:

        if args.taxonomy_spec:

            args.mysql_specs.append(f'taxon={args.taxonomy_spec}')

        if args.all:

            mysql_command = gcm.construct_sql_output_query(args.database_table,
                                                           args.table_columns,
                                                           args.mysql_specs)

        else:

            #Fetch names/db_ids
            query_names = gcm.construct_sql_output_query(None,
                                                         ['name', 'db_id'],
                                                         args.mysql_specs)

            names_dict = gcm.fetch_names(query_names)

            #Fetch primary keys for current versions
            current_ids = gcm.fetch_current_ids(names_dict)

            #Check whether any bioentry tables are being queried, as this
            # affects our new specification
            tables, _, _, _ = gcm.sql_cols(args.database_table,
                                           args.table_columns, args.mysql_specs)

            BIO_TABLES = ['bioentry', 'bioentry_dbxref',
                          'bioentry_qualifier_value',
                          'bioentry_reference', 'biosequence',
                          'comment', 'seqfeature']

            bios = list(set(tables) & set(BIO_TABLES))

            #Construct new specification
            if bios:
                new_spec = [f"(bioentry_id, metadata_id) IN "
                            f"{tuple(current_ids.values())}"]
            else:
                new_spec = [f"metadata_id IN "
                            f"{tuple([ids[1] for ids in current_ids.values()])}"]

            #Construct new SQL query
            mysql_command = gcm.construct_sql_output_query(args.database_table,
                                                           args.table_columns,
                                                           new_spec)

    gcm.csv_from_sql(mysql_command, args.output_name)

elif args.output_format == 'COUNT':

    if args.mysql_query:

        mysql_command = re.sub('SELECT.*?FROM', 'SELECT COUNT(*) FROM',
                               args.mysql_query, 1)

    else:

        if args.taxonomy_spec:

            args.mysql_specs.append(f'taxon={args.taxonomy_spec}')

        if args.all:

            mysql_command = gcm.construct_sql_output_query(None, ['count'],
                                                           args.mysql_specs)

            gcm.return_count(mysql_command)

        else:

            query_names = gcm.construct_sql_output_query(None,
                                                         ['name', 'db_id'],
                                                         args.mysql_specs)

            names_dict = gcm.fetch_names(query_names)

            print(len(names_dict))

else:

    if args.mysql_query:

        mysql_command = args.mysql_query

    else:

        if args.taxonomy_spec:

            args.mysql_specs.append(f'taxon={args.taxonomy_spec}')

        mysql_command = gcm.construct_sql_output_query(None, ['name', 'db_id'],
                                                       args.mysql_specs)

    names_dict = gcm.fetch_names(mysql_command)

    records = gcm.fetch_recs(names_dict, args.all)

    if args.genes:

        records = gcm.extract_genes(records, args.genes)

    gcm.seqfile_from_sql(records, args.output_name, args.output_format.lower())

print('Done.')

print(query_names)
"""
TO-DO

1. Give option to pull either all or just the latest version of each record satisfying spec
    -> If latest, it must run through the master table at some point every time
    -> If latest (i.e. if not --all), consult master table, pull primary keys, load into spec: WHERE 
2. Include ' IN ' in specs?
3. Make sure every query goes through master table


(metadata.metadata_id IN (SELECT master.metadata_id FROM master join metadata on 
master.metadata_id=metadata.metadata_id where <spec>))

((bioentry.bioentry_id, metadata.metadata_id) IN (SELECT master.bioentry_id, master.metadata_id FROM master JOIN metadata ON 
master.metadata_id=metadata.metadata_id JOIN bioentry ON master.bioentry_id=bioentry.bioentry_ID WHERE <spec>))

A) retrieve names satisfying spec       A) Add additional spec
B) retrieve keys from masdter table     B)
C) replace spec                         C)

USAGE NOTES:
1. Assumes taxon name provided is the scientific name
2. GENES: can it be made such that if flag is not given, full genome assumed, if it is given with no argument, all 13 genes are assumes, and if genes specified, then only those genes assumed?
3. Any command entered after -q will be parsed. So be cautious.
4. Outputs most recent by default
5. Assumes user knows table cols
"""

