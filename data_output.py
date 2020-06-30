#!/usr/bin/env python

"""Data retreival from MySQL database.
"""
import argparse
import pandas as pd
import sys

# Import the module:
import gb_csv_module as gcm


#Create the top-level parser
parser = argparse.ArgumentParser(description="Extracting data from the database.")
parser.add_argument('-q', '--query', dest='mysql_query', metavar='{custom_MySQL_query}', help = "Custom MySQL query to extract data from database. (e.g. \"SELECT * FROM metadata WHERE country='United Kingdom';\") NOTE: Your custom specification must be enclosed by double quotations, as above.")
subparsers = parser.add_subparsers(dest="output_format", description='Desired output format:')

#Create the parser for the 'COUNT' command
parser_count = subparsers.add_parser('COUNT', help='Prints an integer on the command line.')
parser_count.add_argument('-s', '--specifications', dest='mysql_specs', metavar='{specifications}', nargs='+', help="Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') REMEMBER: Each individual spec must be enclosed by quotations and separated from the next by a space as above.")

#Create the parser for the 'CSV' command
parser_csv = subparsers.add_parser('CSV', help='Ouputs a .csv file.')
parser_csv.add_argument('-o', '--out', dest = 'output_name', metavar='{Output filename}', required=True, help="Preferred filename for the output (extension will be added automatically according to your output format choice).")
parser_csv.add_argument('-t', '--table', dest='database_table', metavar='{Table name}', choices=["metadata", "bioentry", "bioentry_dbxref", "bioentry_qualifier_value", "bioentry_reference", "biosequence", "seqfeature", "comment", "taxon", "taxon_name"], help="Name of database table you wish to extract data from. (e.g. metadata.) If provided, the script will assume you require data from every column in the specified table. If this is not the case then stating only the required columns under the -c flag will suffice.")
parser_csv.add_argument('-c', '--columns', dest='table_columns', metavar='{Column name(s)}', nargs='+', default=['*'], help="Name of table columns you wish to extract data from. (e.g. name length description)")
parser_csv.add_argument('-s', '--specifications', dest='mysql_specs', metavar='{Specification(s)}', nargs='+', help="Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') NOTE: Each individual specification must be enclosed by quotations and separated from the next by a space as above.")

#Create the parser for the 'FASTA' command
parser_fasta = subparsers.add_parser('FASTA', help='Ouputs a .fasta file.')
parser_fasta.add_argument('-o', '--out', dest = 'output_name', metavar='{Output filename}', required=True, help="Preferred filename for the output (extension will be added automatically according to your output format choice).")
parser_fasta.add_argument('-s', '--specifications', dest='mysql_specs', metavar='{Specification(s)}', nargs='+', help="Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') NOTE: Each individual specification must be enclosed by quotations and separated from the next by a space as above.")

#Create the parser for the 'GB' command
parser_gb = subparsers.add_parser('GB', help="Outputs an annotated .gb file.")
parser_gb.add_argument('-o', '--out', dest = 'output_name', metavar='{Output filename}', required=True, help="Preferred filename for the output (extension will be added automatically according to your output format choice).")
parser_gb.add_argument('-s', '--specifications', dest='mysql_specs', metavar='{Specification(s)}', nargs='+', help="Comma-separated list of mysql specifications.\n(e.g. 'subregion=Sabah' 'length>25000' 'order=Coleoptera') NOTE: Each individual specification must be enclosed by quotations and separated from the next by a space as above.")

# args = parser.parse_args(["-sqlu", "root", "-sqlpw", "mmgdatabase", "-db", "mmg_test", "-t", "metadata", "-o", "metadateru", "-s", "subregion='Sabah'"])
# args = parser.parse_args(["-sqlu", "root", "-sqlpw", "mmgdatabase", "-db", "mmg_test", "-t", "metadata", "-c", "name", "length", "accession", "seq", "-o", "metadateru", "-s", "country='United Kingdon' description='Lucanus sp. BMNH 1425267 mitochondrion, complete genome'", "-f", "csv"])

args=parser.parse_args()

#gcm.csv_from_sql(args.sql_user, args.sql_password, args.db_name, args.database_table, args.table_columns, args.mysql_specs, args.output_name)



if args.output_format == 'CSV':

    if args.mysql_query is None:

        mysql_command = gcm.construct_sql_command(args.database_table, args.table_columns, args.mysql_specs)

    else:

        mysql_command = args.mysql_query

    #gcm.csv_from_sql(mysql_command, args.output_name)

elif args.output_format == 'COUNT':

    if args.mysql_query is None:

        mysql_command = gcm.construct_sql_command(None, ['count'], args.mysql_specs)

    else:

        mysql_command = args.mysql_query

    #gcm.return_count(mysql_command)

else:

    if args.mysql_query is None:

        mysql_command = gcm.construct_sql_command(None, ['name', 'db_id'], args.mysql_specs)

    else:

        mysql_command = args.mysql_query

    """
    names_dict = gcm.fetch_names(mysql_command)

    records = gcm.fetch_recs(names_dict)

    gcm.seqfile_from_sql(records, args.output_name, args.subcommand.lower())
    """


print(f"QUERY: {mysql_command}")


"""
TO-DO

3. What happens if they have chosen multiple columns but have chosen an output format other than CSV?
6. Sort out help page.
7. Number of file output formats needs to be looked at...
9. <>= the only symbols in a query? NO: '!='
10. Data is stored only at the most specific taxonomic level. E.g. a Mordellidae (family) is still a Coleoptera (order), but will only come up under the spec family='Mordellidae' (and not under order='Coleoptera'). 
11. User must input db login details! Not automatic 
12. Add new fields to metadata 

USAGE NOTES:
1. Assumes taxon name provided is the scientific name
2. 

"""
