#!/usr/bin/env python

"""Data retreival from MySQL database. Outputs data in specified format
    according to SQL specifications."""

## Imports ##
import argparse
import sys
import re
import gb_csv_module as gcm
from argparse_formatter import FlexiFormatter

## Arguments ##
parser = argparse.ArgumentParser(
    description="Extracting data from the database.",
    epilog="""
        This script provides 4 output format options, each corresponding \
        to a different subparser: COUNT, CSV, FASTA, GB.
        
        Each subparser contains its own set of arguments and help file. (To \
        see these, simply add the ‘-h’ flag immediately after the parser you \
        want information on — e.g. ‘python3 data_output.py CSV -h’.)
        
        Data is selected from the database according to specifications \
        provided by the user on the command line. Each specification should \
        consist of a field, followed by a comparison operator, followed by \
        the value/parameters the user wishes the output records to satisfy in \
        that field. (e.g. 'subregion=Sabah', 'length>25000', \
        'order=Coleoptera').
        
        For this to be viable, some knowledge of the fields in the database \
        tables is obviously required. The following fields and operators \
        are recognised by this script and can therefore be used in your \
        specifications: 
        
        --FIELDS--
          metadata:\t'metadata_id', 'contigname', 'db_id', 'institution_code', 
                    'collection_code', 'specimen_id', 'morphospecies', \
        'ncbi_taxon_id', 'custom_lineage', 'traptype', 'dev_stage', 'site', \
        'locality', 'subregion', 'country', 'latitude', 'longitude', 'size', \
        'feeding_behaviour', 'habitat', 'habitat_stratum', 'metadata.authors', \
        'library', 'datasubmitter', 'projectname', 'genbank_accession', \
        'notes', 'metadata.version'
        
          bioentry:\t'bioentry_id', 'bioentry.biodatabase_id', 'taxon_id',
                    'name', 'accession', 'identifier', 'division', \
                    'description', 'version'
        
          biosequence:\t'biosequence.version', 'length', 'alphabet', 'seq'
        
          comment:\t'comment_id', 'comment.bioentry_id', 'comment_text',
                   'comment.rank'
                 
          taxon:\t'taxon.taxon_id', 'taxon.ncbi_taxon_id', 'parent_taxon_id',
                 'node_rank', 'genetic_code', 'mito_genetic_code', \
                 'left_value', 'right_value'
               
          taxon_name: 'taxon_name.taxon_id', 'taxon_name.name', 'name_class'
        
        
        --OPERATORS-- 
          =  : Equal to
          != : Not equal to
          <  : Less than
          >  : Greater than
          <= : Less than or equal to
          >= : Greater than or equal to
          IN : TRUE if the operand is equal to one of a list of expressions

        Once parsed, your specifications will be converted into a MySQL \
        'SELECT' query to be executed in the database, returning the result \
        set in your chosen format.
        
        --NOTES-- 
        * For flags accepting multiple arguments (e.g. -s), each argument \
        must be a string in itself. The parser will read spaces as delimiters, \
        so if one of your args contains a space within it (e.g. \
        ‘country=United Kingdom’), it is important that you enclose it with \
        quotation marks to notify the parser that the space is part of one \
        argument string rather than a delimiter between two. Similarly, if \
        your argument string has any single quotation marks in it (e.g. \
        “country IN (‘United Kingdom’, ‘Honduras’)”) , it is important to \
        enclose it with double quotation marks to notify the parser that the \
        single quotes are part of the string rather than delimiters. 
        """,
    formatter_class=FlexiFormatter)
req_group = parser.add_argument_group('required arguments')
req_group.add_argument('--db_user', dest='db_user', help="Database username",
                       required=True, metavar='username')
req_group.add_argument('--db_pass', dest='db_pass', help="Database password",
                       required=True, metavar='password')
parser.add_argument('--all', help="""Use this flag if you wish to pull all \
                    versions of each record satisfying your query. By default \
                    this script will only pull current versions if this flag \
                    is omitted.""", action='store_true')
parser.add_argument('-q',
                    help="""
                    Custom MySQL 'SELECT' query to extract data from \
                    database. 
                    
                    NOTE:
                    1. This must be a correctly-syntaxed MySQL 'SELECT' \
                    statement.
                    2. Your custom specification must be enclosed by double \
                    quotations if it contains any single quotations.
                    
                    EXAMPLE:
                    \"SELECT * FROM metadata WHERE country='United Kingdom';\"
                    """,
                    dest='custom_query')

base_subparser = argparse.ArgumentParser(add_help=False)
base_subparser.add_argument('-s',
                            help="""
                            Space-delimited list of MySQL \
                            specifications for output records to satisfy — 
                            '<FIELD><OPERATOR><VALUE(S)>'. For more \
                            information on hwo to construct these, refer to \
                            the main help page.
                            
                            EXAMPLES:
                              'subregion=Sabah', 
                              'length>25000', 
                              'order=Coleoptera'
                            """,
                            dest='mysql_specs', default=[], nargs='+',
                            metavar='specs')
base_subparser.add_argument('-x',
                            help="""Taxonomic searchterm to search in \
                            database""", dest='taxonomy_spec',
                            metavar='taxon')

subparsers = parser.add_subparsers(dest="output_format",
                                   description='Desired output format:')

#Create the parser for the 'COUNT' command
parser_count = subparsers.add_parser('COUNT', help="""Prints an integer on the \
                                    command line —— For help file see \
                                    'data_output.py COUNT -h'.""",
                                     parents=[base_subparser])

#Create the parser for the 'CSV' command
parser_csv = subparsers.add_parser('CSV', help="""Ouputs a .csv file —— For \
                                    help file see 'data_output.py CSV -h'.""",
                                   parents=[base_subparser])
parser_csv.add_argument('-t', help="""Name of database table you wish
                        to extract data from. (e.g. metadata.) If provided, 
                        the script will assume you require data from every 
                        column in the specified table. If this is not the case 
                        then stating only the required columns under the -c flag
                         will suffice.""", dest='database_table',
                        metavar='table',
                        choices=["metadata", "bioentry", "biosequence",
                                 "taxon", "taxon_name", "rejected", "master"])
parser_csv.add_argument('-c', help="""Name of table columns you 
                        wish to extract data from. (e.g. name length 
                        description)""", dest='table_columns', default=['*'],
                        nargs='+', metavar='columns')
parser_csv.add_argument('-o', help="""Preferred path/filename for the output
                        (extension will be added automatically according to 
                        your output format choice).""", dest='output_name',
                        required=True, metavar='output')

#Create the parser for the 'FASTA' command
parser_fasta = subparsers.add_parser('FASTA', help="""Ouputs a .fasta file —— \
                                    For help file see 'data_output.py FASTA \
                                    -h'.""", parents=[base_subparser])
parser_fasta.add_argument('-g', help="""Name of mitochondrial genes 
                            you wish to extract.""", dest='genes',
                          metavar='genes', nargs='+',
                          choices=['*', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3',
                                   'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L',
                                   'ND5', 'ND6'])
parser_fasta.add_argument('-o', help="""Preferred path/filename for the
                            output (extension will be added automatically 
                            according to your output format choice).""",
                          dest='output_name', metavar='output',
                          required=True)

#Create the parser for the 'GB' command
parser_gb = subparsers.add_parser('GB', help="""Outputs an annotated .gb file \
                                —— For help file see 'data_output.py GB -h'.""",
                                  parents=[base_subparser])
parser_gb.add_argument('-g', help="""Name of mitochondrial genes you 
                        wish to extract.""", dest='genes',
                       metavar='genes', nargs='+',
                       choices=['*', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3',
                                'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L',
                                'ND5', 'ND6'])
parser_gb.add_argument('-o', help="""Preferred path/filename for the 
                        output (extension will be added automatically according 
                        to your output format choice).""", dest='output_name',
                       metavar='output', required=True)

#args = parser.parse_args("--db_user root --db_pass mmgdatabase --all CSV -c order -s 'order=Coleoptera' -o jhgj".split(' '))
#args = parser.parse_args("--db_user root --db_pass mmgdatabase CSV -c db_id ncbi_taxon_id species genus tribe subfamily family superfamily infraorder order class -s country=Honduras -o output/Taxonomy_Testing".split(' '))

args = parser.parse_args()

## Functions ##

# Check login details
gcm.check_login_details(args.db_user, args.db_pass)

# Check any custom queries are 'SELECT' statements
if args.custom_query:
    if not args.custom_query.startswith('SELECT'):
        sys.exit("ERROR: your custom query is not a MySQL 'SELECT' statement.")

## Functions for CSV subparser ##

if args.output_format == 'CSV':

    if args.custom_query:

        if args.taxonomy_spec or args.mysql_specs or args.database_table or \
                args.table_columns:
            # Define restrictions
            parser.error("""query-building specifications and custom mysql \
                query both provided. Please proceed with one method only.""")

        mysql_command = args.custom_query

        gcm.csv_from_sql(mysql_command, args.output_name)

    else:

        # Create list of requested taxonomic fields
        taxreqs = set()

        if args.database_table == 'metadata':
            taxreqs.update(gcm.taxlevels())

        if args.table_columns:
            for col in args.table_columns:
                if col in gcm.taxlevels():
                    taxreqs.add(col)

        if args.taxonomy_spec:
            # Add taxonomic searchterm to specifications if parsed
            args.mysql_specs.append(f'taxon_searchterm={args.taxonomy_spec}')

        if args.table_columns != ['*'] and 'db_id' not in args.table_columns:
            # Add db_id as index for all CSV outputs
            args.table_columns.append('db_id')

        if args.all:

            if taxreqs:

                query_names = gcm.construct_sql_output_query(None,
                                                             ['contigname',
                                                              'db_id'],
                                                             args.mysql_specs)

                names_dict = gcm.fetch_names(query_names)

                current_ids = gcm.fetch_current_ids(names_dict)

            sql = gcm.construct_sql_output_query(args.database_table,
                                                 args.table_columns,
                                                 args.mysql_specs)

            mysql_command = re.sub('SELECT', 'SELECT DISTINCT', sql, 1)

        else:

            # Fetch names/db_ids
            query_names = gcm.construct_sql_output_query(None,
                                                         ['contigname',
                                                          'db_id'],
                                                         args.mysql_specs)

            names_dict = gcm.fetch_names(query_names)

            # Fetch primary keys for current versions
            current_ids = gcm.fetch_current_ids(names_dict)

            # Check whether any bioentry tables are being queried, as this
            # affects our new specification
            tables, _, _, _ = gcm.sql_cols(args.database_table,
                                           args.table_columns, args.mysql_specs)

            BIO_TABLES = ['bioentry', 'bioentry_dbxref',
                          'bioentry_qualifier_value',
                          'bioentry_reference', 'biosequence',
                          'comment', 'seqfeature']

            bios = list(set(tables) & set(BIO_TABLES))

            # Construct new specification
            if bios:
                if len(current_ids) == 1:
                    # Spec must be built without tuple as single-element
                    # tuples contain a trailing comma which breaks sql query.
                    new_spec = [f"(bioentry_id, metadata_id)="
                                f"{list(current_ids.values())[0]}"]
                else:
                    new_spec = [f"(bioentry_id, metadata_id) IN " 
                                f"{tuple(current_ids.values())}"]
            else:
                if len(current_ids) == 1:
                    new_spec = [f"metadata_id="
                                f"{list(current_ids.values())[0][1]}"]
                else:
                    meta_ids = tuple([ids[1] for ids in current_ids.values()])
                    new_spec = [f"metadata_id IN {meta_ids}"]

            # Construct new SQL query
            mysql_command = gcm.construct_sql_output_query(args.database_table,
                                                           args.table_columns,
                                                           new_spec)

        df_out = gcm.df_from_sql(mysql_command)

        if taxreqs:
            # Create taxonomy dict
            taxonomy = {dbid: gcm.fetch_taxonomy(current_ids[dbid][0], taxreqs)
                        for dbid in current_ids.keys()}

            # Add taxonomy cols to dataframe
            # TODO: Reorder dataframe columns to put them in x order in x part of the header
            gcm.add_taxonomy_to_df(df_out, taxonomy, taxreqs)

        # Delete surrogate keys
        df_out.drop(['version', 'metadata_id', 'bioentry_id'], errors='ignore',
                    inplace=True)

        # Write file
        gcm.csv_from_df(df_out, args.output_name)

    print('Done.')

## Functions for COUNT subparser ##

elif args.output_format == 'COUNT':

    if args.custom_query:

        if args.taxonomy_spec or args.mysql_specs:
            # Define restrictions
            parser.error("""query-building specifications and custom mysql \
                query both provided. Please proceed with one method only.""")

        if args.all:

            #One would have to do a normal SELECT command here
            mysql_command = re.sub('SELECT.*?FROM',
                                   'SELECT DISTINCT COUNT(*) FROM',
                                   args.custom_query, 1)

            gcm.return_count(mysql_command)

        else:

            mysql_command = re.sub('SELECT.*?FROM',
                                   'SELECT metadata.contigname, '
                                   'metadata.db_id FROM',
                                   args.custom_query, 1)

            names = gcm.fetch_names(mysql_command)

            print(len(names))

    else:

        if args.taxonomy_spec:

            args.mysql_specs.append(f'taxon_searchterm={args.taxonomy_spec}')

        if args.all:

            sql = gcm.construct_sql_output_query(None, ['count'],
                                                 args.mysql_specs)

            mysql_command = re.sub('SELECT', 'SELECT DISTINCT', sql, 1)

            gcm.return_count(mysql_command)

        else:

            query_names = gcm.construct_sql_output_query(None,
                                                         ['contigname',
                                                          'db_id'],
                                                         args.mysql_specs)

            names = gcm.fetch_names(query_names)

            print(len(names))

## Functions for FASTA and GB subparsers ##

else:

    if args.custom_query:

        if args.taxonomy_spec or args.mysql_specs:
            # Define restrictions
            parser.error("""query-building specifications and custom mysql \
                query both provided. Please proceed with one method only.""")

        mysql_command = args.custom_query

    else:

        if args.taxonomy_spec:

            args.mysql_specs.append(f'taxon_searchterm={args.taxonomy_spec}')

        mysql_command = gcm.construct_sql_output_query(None, ['contigname',
                                                              'db_id'],
                                                       args.mysql_specs)

    names_dict = gcm.fetch_names(mysql_command)

    records = gcm.fetch_recs(names_dict, args.all)

    if args.genes:

        records = gcm.extract_genes(records, args.genes)

    gcm.seqfile_from_sql(records, args.output_name, args.output_format.lower())

    print('Done.')