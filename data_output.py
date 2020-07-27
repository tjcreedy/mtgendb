#!/usr/bin/env python

"""Data retreival from MySQL database.
"""
import argparse
import re

#Import the module:
import gb_csv_module as gcm


#Create the top-level parser
parser = argparse.ArgumentParser(description="Extracting data from the database.")
req_group = parser.add_argument_group('Required arguments')
req_group.add_argument('--db_user', dest='db_user', help="Database username", metavar='{db_username}', required=True)
req_group.add_argument('--db_pass', dest='db_pass', help="Database password", metavar='{db_password}', required=True)
parser.add_argument('-q', '--query', dest='mysql_query', help="""Custom MySQL 
query to extract data from database. (e.g. \"SELECT * FROM metadata WHERE 
country='United Kingdom';\") NOTE: Your custom specification must be enclosed by double quotations, as above.""",
                    metavar='custom_MySQL_query')
subparsers = parser.add_subparsers(dest="output_format", description='Desired output format:')

#Create the parser for the 'COUNT' command
parser_count = subparsers.add_parser('COUNT', help="""Prints an integer on the command line —— 
For help file see 'data_output.py COUNT -h'.""")
parser_count.add_argument('-s', '--specifications', dest='mysql_specs', help="""Comma-separated 
list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') 
REMEMBER: Each individual spec must be enclosed by quotations and separated from the 
next by a space as above.""", metavar='{specifications}', nargs='+')

#Create the parser for the 'CSV' command
parser_csv = subparsers.add_parser('CSV', help="""Ouputs a .csv file —— For 
help file see 'data_output.py CSV -h'.""")
parser_csv.add_argument('-o', '--out', dest = 'output_name',
                        help="""Preferred filename for the output (extension 
                        will be added automatically according to your 
                        output format choice).""",
                        metavar='{Output filename}', required=True)
parser_csv.add_argument('-t', '--table', dest='database_table', metavar='{Table name}', choices=["metadata", "bioentry", "bioentry_dbxref", "bioentry_qualifier_value", "bioentry_reference", "biosequence", "seqfeature", "comment", "taxon", "taxon_name"], help="Name of database table you wish to extract data from. (e.g. metadata.) If provided, the script will assume you require data from every column in the specified table. If this is not the case then stating only the required columns under the -c flag will suffice.")
parser_csv.add_argument('-c', '--columns', dest='table_columns', metavar='{Column name(s)}', nargs='+', default=['*'], help="Name of table columns you wish to extract data from. (e.g. name length description)")
parser_csv.add_argument('-s', '--specifications', dest='mysql_specs', metavar='{Specification(s)}', nargs='+', help="Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') NOTE: Each individual specification must be enclosed by quotations and separated from the next by a space as above.")

#Create the parser for the 'FASTA' command
parser_fasta = subparsers.add_parser('FASTA', help="Ouputs a .fasta file —— For help file see 'data_output.py FASTA -h'.")
parser_fasta.add_argument('-o', '--out', dest = 'output_name', metavar='{Output filename}', required=True, help="Preferred filename for the output (extension will be added automatically according to your output format choice).")
parser_fasta.add_argument('-s', '--specifications', dest='mysql_specs', metavar='{Specification(s)}', nargs='+', help="Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') NOTE: Each individual specification must be enclosed by quotations and separated from the next by a space as above.")
parser_fasta.add_argument('-g', '--genes', dest='genes', metavar='{gene_names}', nargs='+', choices=['*', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6'], help='Name of mitochondrial genes you wish to extract.')

#Create the parser for the 'GB' command
parser_gb = subparsers.add_parser('GB', help="Outputs an annotated .gb file —— For help file see 'data_output.py GB -h'.")
parser_gb.add_argument('-o', '--out', dest = 'output_name', metavar='{Output filename}', required=True, help="Preferred filename for the output (extension will be added automatically according to your output format choice).")
parser_gb.add_argument('-s', '--specifications', dest='mysql_specs', metavar='{Specification(s)}', nargs='+', help="Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') NOTE: Each individual specification must be enclosed by quotations and separated from the next by a space as above.")
parser_gb.add_argument('-g', '--genes', dest='genes', metavar='{gene_names}', nargs='+', choices=['*', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6'], help='Name of mitochondrial genes you wish to extract.')

# args = parser.parse_args(["-sqlu", "root", "-sqlpw", "mmgdatabase", "-db", "mmg_test", "-t", "metadata", "-o", "metadateru", "-s", "subregion='Sabah'"])
# args = parser.parse_args(["-sqlu", "root", "-sqlpw", "mmgdatabase", "-db", "mmg_test", "-t", "metadata", "-c", "name", "length", "accession", "seq", "-o", "metadateru", "-s", "country='United Kingdon' description='Lucanus sp. BMNH 1425267 mitochondrion, complete genome'", "-f", "csv"])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'GB', '-o', 'replace', '-s',  'country=China'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'FASTA', '-o', 'outksis', '-s',  'country!=Malaysia', '-g' , 'COX2', 'ND3', 'ATP6'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'FASTA', '-o', 'outksis', '-s',  'country!=Malaysia', '-g' , '*'])
# args = parser.parse_args(['--db_user', 'root', '--db_pass', 'mmgdatabase', 'FASTA', '-o', 'malays', '-s',  'species=Stenus boops', '-g', 'COX2', 'ND3', 'ATP6'])


args = parser.parse_args()

#gcm.csv_from_sql(args.sql_user, args.sql_password, args.db_name, args.database_table, args.table_columns, args.mysql_specs, args.output_name)

if args.output_format == 'CSV':

    if not args.mysql_query:

        mysql_command = gcm.construct_sql_output_query(args.database_table, args.table_columns, args.mysql_specs)

    else:

        mysql_command = args.mysql_query

    gcm.csv_from_sql(mysql_command, args.output_name, args.db_user, args.db_pass)

elif args.output_format == 'COUNT':

    if not args.mysql_query:

        mysql_command = gcm.construct_sql_output_query(None, ['count'], args.mysql_specs)

    else:

        mysql_command = re.sub('SELECT.*?FROM', 'SELECT COUNT(*) FROM', args.mysql_query, 1)

    gcm.return_count(mysql_command, args.db_user, args.db_pass)

else:

    if not args.mysql_query:

        mysql_command = gcm.construct_sql_output_query(None, ['name', 'db_id'], args.mysql_specs)

    else:

        mysql_command = args.mysql_query

    names_dict = gcm.fetch_names(mysql_command, args.db_user, args.db_pass)

    records = gcm.fetch_recs(names_dict, args.db_user, args.db_pass)

    if args.genes:

        records = gcm.extract_genes(records, args.genes)

    gcm.seqfile_from_sql(records, args.output_name, args.output_format.lower())


#print(f"QUERY: {mysql_command}")


"""
TO-DO

7. Number of file output formats... GENBANK SUBMISSION - .sqn?

USAGE NOTES:
1. Assumes taxon name provided is the scientific name
2. GENES: can it be made such that if flag is not given, full genome assumed, if it is given with no argument, all 13 genes are assumes, and if genes specified, then only those genes assumed?
3. Any command entered after -q will be parsed. So be cautious.
"""

