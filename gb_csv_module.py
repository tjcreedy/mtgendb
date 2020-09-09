#!/urs/bin/env python

"""Functions to interact with data files and MySQL database."""


## Imports ##
import sys, time, urllib.request, csv, re
from collections import defaultdict
from Bio import SeqIO, Entrez
from Bio.SeqFeature import SeqFeature, FeatureLocation

import pandas as pd
import numpy as np

import MySQLdb as mdb
from BioSQL import BioSeqDatabase
from sqlalchemy import create_engine, update

# TODO: TEST FOLLOWING COMMANDS:
# 1. Create database, load schema
# 2. gb/csv ingest
# 3. gb ingest
# 4. output chinas. Output chinas and honduras'. output chinas with length > 17000
#    ouptput using own command

## Variables ##
# TODO: Should these begin with a trailing underscore to denote they are for internal use?
db_driver = "MySQLdb"
db_passwd = "mmgdatabase"
db_host = "localhost"
db_user = "root"
db_name = "mmg_test"
mysql_engine = "mysql+mysqldb://root:mmgdatabase@localhost/mmg_test"
namespace = "mmg"

## Functions ##

def check_login_details(db_un, db_pw):
    """Checks database login details are correct.
    """
    if db_un != db_user or db_pw != db_passwd:
        sys.exit("ERROR: Incorrect login details.")

    return

def gb_into_dictionary(gb_filename, key):
    """Take a file in genbank-format and load it into a dictionary, define
    function for keys (default: by name in LOCUS).
    """
    global gb_dictionary

    if key == "LOCUS":

        def get_seqname(record):
            # Given a SeqRecord, return the sequence name as a string.
            seqname = str(record.name)
            return seqname

        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"),
                                      key_function=get_seqname)

    elif key == "ACCESSION":

        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"))

    elif key == "DEFINITION":

        def get_definition(record):
            # Given a SeqRecord, return the definition as a string.
            definition = str(record.description)
            return definition

        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"),
                                      key_function=get_definition)

    return gb_dictionary

def text_to_list(txtfile):
    """Converts text-file of IDs (one per line) into a list, exiting if
    duplicates are present.
    """
    with open(txtfile, "r") as txt:
        # Create a comma-delimited list of accession numbers from text file
        # (stripping any blank spaces/empty lines).
        ids = txt.read()
        striplist = lambda lis: [x.strip() for x in lis]
        ids_list = list(filter(None, striplist(ids.split('\n'))))

    duplicates = set([idd for idd in ids_list if ids_list.count(idd) > 1])

    if len(duplicates):
        # If there are any duplicates, drop these proceed with unque accessions
        # only.
        print("WARNING: Duplicates detected in text file: '" + "', '".join(
            duplicates) + "'. Proceeding with unique entries only.")
        ids_list = list(set(ids_list))

    return ids_list

def versions_to_dict(txtfile):
    """Converts 3-column text-file of IDs and target versions for rollback
    into a dict.

    Argument:
     - 3-column tab-delimited text-file: db_id <TAB> bio_vers <TAB> meta_vers

    Output:
     - Nested dict: {db_id: {'b': bio_vers, 'm': meta_vers}, ... }

    """
    with open(txtfile, 'r') as txt:
        # Create a comma-delimited list of accession numbers from text file
        # (stripping any blank spaces/empty lines).
        txt_string = txt.read()
        lines_raw = list(filter(None, txt_string.split('\n')))
        lines_stripped = [line.strip() for line in lines_raw]

    target_versions = {}
    for line in lines_stripped:
        if line.count('\t') not in [1, 2]:
            sys.exit("ERROR: your input text file must be a 3-column "
                     "tab-delimited list: {db_id} <TAB> {bio_version} <TAB> "
                     "{meta_version}")
        elif line.count('\t') == 1:
            # Only bio_version column provided
            vals = line.split('\t')
            target_versions[vals[0]] = {'b': int(vals[1]), 'm': None}
        else:
            vals = line.split('\t')
            if '' in vals:
                target_versions[vals[0]] = {'b': None, 'm': int(vals[2])}
            else:
                target_versions[vals[0]] = {'b': int(vals[1]),
                                            'm': int(vals[2])}

    return target_versions

def check_acc_format(acc_list):
    """Checks each accession in list satisfies genbank format, dropping any
    that don't.
    """
    accs_out = []

    for acc in acc_list:
        z = re.match("[A-Z]{1,3}_?[0-9]{5,8}$", acc)
        if z:
            accs_out.append(acc)
        else:
            print(f"WARNING: Accession '{acc}' dropped due to incorrect "
                  f"GenBank format.")

    return accs_out

def correct_header(csv_dataframe, action):
    """Check metadata file has correct columns by creating a list of its column
     headers and checking it against a list of the expected headers.

     Arguments:
         - csv_dataframe - dataframe
         - action - ['ingest', 'update', 'rollback', 'remove']
    """
    csv_header = csv_dataframe.columns.values.tolist()

    if action == 'ingest':
        expected_header = ['contigname', 'institution_code', 'collection_code',
                           'specimen_id', 'morphospecies', 'species', 'genus',
                           'subfamily', 'family', 'order', 'taxid', 'traptype',
                           'dev_stage', 'site', 'locality', 'subregion',
                           'country', 'latitude', 'longitude', 'size',
                           'feeding_behaviour', 'habitat', 'habitat_stratum',
                           'authors', 'library', 'datasubmitter', 'projectname',
                           'genbank_accession', 'notes']
    else:
        expected_header = ['contigname', 'db_id', 'institution_code',
                           'collection_code', 'specimen_id', 'morphospecies',
                           'ncbi_taxon_id', 'custom_lineage', 'traptype',
                           'dev_stage', 'site', 'locality', 'subregion',
                           'country', 'latitude', 'longitude', 'size',
                           'feeding_behaviour', 'habitat', 'habitat_stratum',
                           'authors', 'library', 'datasubmitter', 'projectname',
                           'genbank_accession', 'notes']

    if expected_header != csv_header:
        print("Incorrect header in CSV file.\n")
        sys.exit("Current header is: " + str(csv_header) +
                 "\n\nIt must be as follows: " + str(expected_header))

    return

def lat_long(df):
    """Checks the latitude and longitude columns of the DataFrame are correctly
    formatted.
    """
    vals = list(zip(df['contigname'], df['latitude'], df['longitude']))

    for (name, lat, long) in vals:
        if not isinstance(lat, float):
            sys.exit(f"ERROR: latitude for entry '{name}' not of type 'float':"
                     f" '{lat}'")
        if not isinstance(long, float):
            sys.exit(f"ERROR: longitude for entry '{name}' not of type 'float':"
                     f" '{long}'")

    return

def matching_inputids(csv_df, gb_dict, action):
    """Check if the GenBank and CSV metadata files have matching entries, and
    that there are no duplicates in the CSV file — duplicates in the GenBank
    file would have already been flagged in gb_into_dictionary().
    """
    new_gb_dict, new_csv_df = [gb_dict, csv_df]

    if action == 'ingest':
        ids_csv = csv_df['contigname'].values.tolist()
    else:
        ids_csv = csv_df['db_id'].values.tolist()

    union = list(set(ids_csv + list(gb_dict.keys())))
    intersect = list(set(ids_csv) & set(gb_dict.keys()))
    discrepant_ids = [a for a in union if a not in intersect]
    csv_duplicates = set([idd for idd in ids_csv if ids_csv.count(idd) > 1])

    if len(csv_duplicates):
        sys.exit("ERROR: There are multiple rows sharing the same name in your "
                 "CSV file: " + ', '.join(csv_duplicates) + ". IDs must be unique.")

    if len(intersect) != len(union):
        # If the IDS in the CSV and GenBank files are not identical...
        x = input("WARNING: Your CSV and GenBank files contain different "
                  "entries.\nWould you like to ignore these and proceed with "
                  "shared entries only or cancel the operation? "
                  "'P'/'C'\n?>").capitalize()
        while not (x == 'C' or x == 'P'):
            x = input("Type 'P' to ignore discrepant entries and proceed, or "
                      "type 'C' to cancel the operation.\n?>").capitalize()
        if x == 'C':
            sys.exit("Operation cancelled.")
        else:
            # Return new gb_dict with discrepant entries deleted
            for gb_record in gb_dict:
                if gb_record in discrepant_ids:
                    print(f" - Skipping entry '{gb_record}' as it appears in "
                          f"the GenBank file but not the CSV file.")

            new_gb_dict = {key: gb_dict[key] for key in gb_dict if
                           key not in discrepant_ids}

            # Return new csv_df with discrepant entries deleted
            for contigname in ids_csv:
                if contigname in discrepant_ids:
                    print(f" - Skipping entry '{contigname}' as it appears in "
                          f"the CSV file but not the GenBank file.")

            df_missing_ids = [b for b in discrepant_ids if b not in gb_dict]
            new_csv_df.set_index('contigname', inplace=True)
            new_csv_df.drop(df_missing_ids, inplace=True)
            new_csv_df.reset_index(inplace=True)

    return new_csv_df, new_gb_dict

def new_ids(genbank_dict, prefix, startvalue, padding):
    """Check if the new ids are looking fine given the user input (prefix,
    startvalue, padding).
    """
    dict_new_ids = {}

    if len(prefix) > 4:
        sys.exit("ERROR: The prefix (-p) is too long (it should be no more than"
                 " 4 characters).")

    if not prefix.isalpha():
        sys.exit("ERROR: The prefix (-p) should only consist of letters.")

    if not str(startvalue).isdigit():
        sys.exit("ERROR: The number (-n) should only consist of digits.")

    if len(str(len(genbank_dict.keys()) + int(startvalue) - 1)) > int(padding):
        # Checks there's not too many sequences to assign serial numbers to
        # given the padding.
        sys.exit("ERROR: A 0-padding by " + str(padding) +
                 " digits is not sufficient for the ingestion of " + str(
                len(genbank_dict)) + " new entries.")

    if len(str(startvalue)) > int(padding):
        sys.exit("ERROR: The starting number " + str(
            startvalue) + " exceeds the digits for 0-padding (" + str(
            padding) + " digits).")

    for gb_record, record in genbank_dict.items():
        record_name = record.name
        pad = "{0:0" + str(padding) + "d}"
        padded = pad.format(int(startvalue))
        new_name = prefix + str(padded)
        dict_new_ids[record_name] = new_name
        startvalue = int(startvalue) + 1

    return dict_new_ids

def check_ids(ids_list, action):
    """Check if the database ids already exist in the database.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)

    for idd in ids_list:
        mysql_command = f"SELECT * FROM metadata WHERE (contigname='{idd}') " \
                        f"OR (db_id='{idd}');"
        with con:
            cur = con.cursor()
            cur.execute(mysql_command)
            result = cur.fetchone()

        if action in ['replace', 'remove', 'rollback']:
            if result is None:
                sys.exit(f"ERROR: The id '{idd}' is missing from the database.")

        elif action == 'ingest':
            if result is not None:
                sys.exit(f"ERROR: The id '{idd}' already exists in the "
                         f"database:\n{str(result)}")

        else:
            sys.exit(f"ERROR: Unrecognized action: '{action}'")

    return

def check_accs_in_db(accs_list):
    """Checks accessions haven't already been pulled from GenBank by checking
    them against the 'name' column in the metadata table.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)

    duplicates = set()
    for acc in accs_list:
        sql_1 = f"SELECT * FROM metadata WHERE contigname='{acc}';"
        sql_2 = f"SELECT * FROM rejected WHERE accession='{acc}';"
        with con:
            cur = con.cursor()
            cur.execute(sql_1)
            if cur.fetchone():
                duplicates.add(acc)
            cur.execute(sql_2)
            if cur.fetchone():
                duplicates.add(acc)

    if duplicates:
        x = input(f"WARNING: There are accessions in your input file that have"
                  f" previously been rejected or already exist in the database:"
                  f" {', '.join(list(duplicates))}\nWould you like to drop"
                  f" these and proceed or cancel the operation? 'P'/'C'").upper()
        while not (x == 'P' or x == 'C'):
            x = input(f"Type 'P' to drop duplicate accessions and"
                      f" proceed or 'C' to cancel the operation.").upper()
        if x == 'C':
            sys.exit('Operation cancelled.')
        else:
            accs_list = [a for a in accs_list if a not in duplicates]

    return accs_list

def check_recs_for_dups(records):
    """Checks a list of records for duplicates by referring to the comment
    section and comparing seqs.
    """
    recs = list(records.values())
    dupgroups = dict()
    nrecs = len(recs)
    grpn = 0

    for i, seqr in enumerate(recs):
        if seqr.name in dupgroups:
            continue
        dup = {seqr.name}
        if seqr.name.startswith('NC_'):
            # Checks annotation section to see if it's copy of anything.
            # If it is, then it adds its counterpart to dup
            origreg = 'The reference sequence is identical to ([A-Z]+[0-9]+)\.'
            origsearch = re.search(origreg, seqr.annotations['comment'])
            if origsearch:
                # .group(1) grabs parentheses-enclosed regex in above string
                dup.add(origsearch.group(1))
        if i < nrecs:
            for comprec in recs[i + 1:]:
                # Search remaining recs for names in dups, and check
                # for duplicate seqs
                if comprec.name in dup:
                    continue
                if str(seqr.seq) == str(comprec.seq):
                    dup.add(comprec.name)
        for n in dup:
            dupgroups[n] = grpn
        grpn += 1

    groupdups = defaultdict(list)
    for name, group in dupgroups.items():
        groupdups[group].append(name)

    use, reject = [[], []]
    for names in groupdups.values():
        if len(names) == 1:
            use.append(names[0])
        else:
            # If any begin with 'NC_' then use the first of those, otherwise
            # use any.
            nc_i = [n for n in names if n.startswith('NC_')]
            if len(nc_i) > 0:
                use.append(nc_i[0])
                rej = [n for n in names if n != nc_i[0]]
            else:
                use.append(names[0])
                rej = names[1:]

            reject.extend(rej)

    if reject:
        print(f"WARNING: The following records were found to be duplicates of "
              f"other records have been dropped from your"
              f" GenBank records:\n{', '.join(reject)}")

    use_recs = {acc: rec for acc, rec in records.items() if acc in use}

    return use_recs, reject

def check_seqs_in_db(records):
    """Checks whether sequences already exist in the database.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    reject = []
    for acc, rec in records.items():
        sql = f"SELECT * FROM biosequence WHERE seq='{str(rec.seq)}';"
        with con:
            cur = con.cursor()
            cur.execute(sql)
            if cur.fetchone():
                reject.append(acc)

    if reject:
        print(f"WARNING: The following sequences have been dropped as they "
              f"already exist in the database:"
              f"\n{', '.join(reject)}")

    new_recs = {name: rec for name, rec in records.items() if name not in reject}

    return new_recs, reject

def load_reject_table(rejected):
    """Loads a list of rejected accession numbers into the 'rejected' SQL
    table.
    """
    rej_df = pd.DataFrame({'accession': rejected})
    engine = create_engine(mysql_engine, echo=False)
    rej_df.to_sql(name='rejected', if_exists='append', index=False, con=engine)

    return

def check_latest_version(db_id):
    """Return most recent bioentry and metaentry version numbers of
    record
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)

    with con:
        cur = con.cursor()
        sql = f"""SELECT MAX(bioentry.version) , MAX(metadata.version) FROM 
            metadata JOIN bioentry ON metadata.db_id=bioentry.name WHERE 
            metadata.db_id='{db_id}';"""
        cur.execute(sql)
        for row in cur.fetchall():
            bioentry_version, metadata_version = row

    return bioentry_version, metadata_version

def check_current_version(db_id):
    """Return current bioentry and metaentry version numbers of record
    (as represented in the master).
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)

    with con:
        cur = con.cursor()
        sql_ids = f"""SELECT bioentry_id, metadata_id FROM master WHERE 
                db_id='{db_id}';"""
        cur.execute(sql_ids)
        for row in cur.fetchall():
            bioentry_id, metadata_id = row
        sql_versions = f"""SELECT bioentry.version, metadata.version FROM 
            metadata JOIN bioentry ON metadata.db_id=bioentry.name WHERE 
            metadata.metadata_id={metadata_id} AND bioentry.bioentry_id=
            {bioentry_id};"""
        cur.execute(sql_versions)
        for row in cur.fetchall():
            bio_version, meta_version = row

    return bio_version, meta_version

def fetch_current_ids(names_dict):
    """Fetch primary keys of current versions from master table.

    Argument:
     - names_dict - dict with names as keys and db_ids as values
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd,
                      db=db_name)
    current_ids = {}

    with con:
        cur = con.cursor()
        for db_id in names_dict.values():
            sql_ids = f"""SELECT bioentry_id, metadata_id FROM master WHERE 
                db_id='{db_id}';"""
            cur.execute(sql_ids)
            for row in cur.fetchall():
                bioentry_id, metadata_id = row
            current_ids[db_id] = (bioentry_id, metadata_id)

    return current_ids

def fetch_names(sql, db_un, db_pw):
    """Fetch names and corresponding db_id's from database using MySQL command

    Argument:
     - sql - SQL query: "SELECT metadata.contigname, metadata.db_id FROM... ;"
    """
    con = mdb.connect(host=db_host, user=db_un, passwd=db_pw, db=db_name)

    with con:
        cur = con.cursor()
        cur.execute(sql)
        records = cur.fetchall()

    names_dict = {row[0]: row[1] for row in set(records)}

    return names_dict

def change_ids_genbank(genbank_dict, dict_new_ids, key):
    """Change the ids in the genbank file to the new database ids (LLLLNNNNN).
    """
    for gb_record, record in genbank_dict.items():

        record.name = dict_new_ids[record.name]
        new_accession = record.name + ".0"
        record.id = new_accession
        record.annotations["accessions"] = [record.name]
        record.annotations["sequence_version"] = 0
        if key == "DESCRIPTION":
            # Erase the description as it contained the input id
            record.description = ""

    return

def change_names_csv(csv_df, dict_new_ids):
    """Adds column containing new db_ids to df.
    """
    dict_df = pd.DataFrame(list(dict_new_ids.items()),
                           columns=["contigname", "db_id"])
    new_csv_df = pd.merge(dict_df, csv_df, on='contigname')

    return new_csv_df

def taxonomy_metadata(csv_df):
    """Function returning a dictionary with all the taxonomic information for
    each id from metadata ([] if no info given).
    """
    csv_taxa = {}
    ids_csv = csv_df['contigname'].values.tolist()
    df = csv_df.set_index("contigname")

    for id_ in ids_csv:

        tax_list = []
        species = df.loc[id_, 'species']
        subfamily = df.loc[id_, 'subfamily']
        family = df.loc[id_, 'family']
        order = df.loc[id_, 'order']

        if not pd.isnull(species):
            tax_list.append(species)
        if not pd.isnull(subfamily):
            tax_list.append(subfamily)
        if not pd.isnull(family):
            tax_list.append(family)
        if not pd.isnull(order):
            tax_list.append(order)

        csv_taxa[id_] = tax_list

    return csv_taxa

def chunker(seq, size):
    """Splits input list/set into subsets of specified size
    """
    if type(seq) is set:
        seq = list(seq)

    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def return_gb_data(acc_list, email):
    """Takes list of of IDs/accessions and returns dict of corresponding
    GenBank entries.
    """
    print("Fetching GenBank records...")

    records = {}
    Entrez.email = email

    for accset in chunker(acc_list, 200):
        # Fetch gb entries in batches of 200
        with Entrez.efetch(db='nuccore', id=accset, rettype='gb') as handle:
            for record in SeqIO.parse(handle, "gb"):
                records[record.name] = record
        time.sleep(0.1)

    for acc in acc_list:
        if acc not in records.keys():
            print(f" - Accession '{acc}' returned no hits on NCBI.")

    return records

def extract_metadata(records):
    """Extracts metadata from gb_dict and writes to DataFrame.
    """
    # Extract metadata to separate dict
    gb_metadata = {}
    for acc, record in records.items():
        for feature in record.features:
            if feature.type == "source":
                record_metadata = []

                if 'db_xref' in feature.qualifiers:
                    # Extract taxon_id
                    taxid = feature.qualifiers['db_xref'][0].split(":")[1]
                else:
                    print(f"WARNING: 'db_xref' absent from '{acc}' qualifiers.")
                    taxid = None

                if 'country' in feature.qualifiers:
                    # Extract country
                    country = feature.qualifiers['country'][0]
                    if ': ' in country:
                        country = feature.qualifiers['country'][0].split(": ")[0]
                        subregion = feature.qualifiers['country'][0].split(": ")[1]
                    else:
                        subregion = None
                else:
                    country, subregion = [None, None]

                if 'isolation_source' in feature.qualifiers:
                    # Extract collection_method
                    traptype = feature.qualifiers['isolation_source'][0]
                else:
                    traptype = None

                if 'specimen_voucher' in feature.qualifiers:
                    # Extract accession, stripping all instances of ' ', '-',
                    # '.' and ':' from string.
                    accession = feature.qualifiers['specimen_voucher'][
                        0].translate({ord(c): None for c in
                                      ' -.:'})
                else:
                    accession = None

                record_metadata.append(taxid)
                record_metadata.append(country)
                record_metadata.append(subregion)
                record_metadata.append(traptype)
                record_metadata.append(accession)

                gb_metadata[acc] = record_metadata

    # Write to DataFrame
    gb_met_df = pd.DataFrame.from_dict(gb_metadata, orient='index')
    gb_met_df.reset_index(inplace=True)
    gb_met_df.columns = ["contigname", "ncbi_taxon_id", "country", "subregion",
                         "traptype", "genbank_accession"]

    # Add other metadata columns required by database and fill them with blanks
    # (None/NaN objects)
    for label in ['institution_code', 'collection_code', 'specimen_id',
                  'morphospecies', 'custom_lineage', 'dev_stage', 'site',
                  'locality', 'authors', 'size', 'habitat', 'feeding_behaviour',
                  'locomotion', 'library', 'datasubmitter',
                  'projectname', 'notes']:
        gb_met_df[label] = None

    for header in ['latitude', 'longitude']:
        gb_met_df[header] = np.nan

    return gb_met_df

def taxonomy_from_gb(genbank_dict):
    """Function returning dictionary of taxids from genbank file if it has it
    (or "" if not).
    """
    gb_taxa = {}

    for gb_record, record in genbank_dict.items():
        tax_id = ""
        for (index, feature) in enumerate(record.features):
            if feature.type == "source":
                if "db_xref" in feature.qualifiers:
                    xref = feature.qualifiers["db_xref"]

                    if len(xref) > 1:
                        print("There are several ids in the db_xref 'source' "
                              "information.")
                    else:
                        tax_id = xref[0].split(":")[1]

        gb_taxa[record.name] = tax_id

    return gb_taxa

def taxid_metadata(csv_dataframe):
    """Return dictionary with taxids provided in metadata table
    ("" if not given).
    """
    csv_taxids = {}
    ids_csv = csv_dataframe['contigname'].values.tolist()
    df = csv_dataframe.set_index("contigname")

    for id_ in ids_csv:  # for each id in ids_csv...

        taxid = df.loc[id_, 'taxid']

        if not pd.isnull(taxid):
            csv_taxids[id_] = taxid
        else:
            csv_taxids[id_] = ""

    return csv_taxids

def return_ncbi_taxid(entry, searchterm, email_address):
    """For each entry get tax_id from NCBI taxonomy based on taxonomic
    information.
    """
    Entrez.email = email_address
    handle = Entrez.esearch(db="taxonomy", retmax=2, term=searchterm)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]

    if len(id_list) == 0:
        # Give user option to use unsuccessful searchterm as custom lineage info
        # or cancel the operation.
        x = input(f" - No hits found for search term '{searchterm}' in NCBI "
                  f"taxonomy.\n   Would you like to record this as custom "
                  f"lineage information for entry '{entry}' and proceed or "
                  f"cancel the operation? 'P'/'C'\n ?>").capitalize()
        while not (x == 'P' or x == 'C'):
            x = input(f"   Type 'P' to record '{searchterm}' as custom lineage "
                      f"information for entry '{entry}' or 'C' to cancel the "
                      f"operation.\n ?>").capitalize()
        if x == 'C':
            sys.exit("\nOperation cancelled.")
        else:
            tax_id = ""
            print(f" - '{searchterm}' saved to custom lineage information for "
                  f"entry '{entry}'.\n...")

    elif len(id_list) > 1:
        print(f" - Multiple hits found for search term '{searchterm}' in NCBI "
              f"taxonomy.")
        tax_id = ""

    elif len(id_list) == 1:
        tax_id = id_list[0]

    time.sleep(0.5)

    return tax_id

def return_ncbi_lineage(searchterm, email_address):
    """Search NCBI for lineage information given a tax id.
    """
    Entrez.email = email_address
    handle = Entrez.efetch(db="taxonomy", id=searchterm)
    record = Entrez.read(handle)
    handle.close()

    if len(record) == 0:
        print(f"No hits found for tax id '{searchterm}' in NCBI taxonomy.")

    elif len(record) > 1:
        print(f"Multiple hits found for tax id '{searchterm}' in NCBI taxonomy.")

    elif len(record) == 1:
        taxonomy = record[0]["Lineage"]
        taxon = record[0]["ScientificName"]
        lineage = taxonomy + "; " + taxon

    time.sleep(0.5)

    return lineage


def get_ncbi_lineage(csv_dataframe, email_address, searchterm):
    """Search on NCBI for tax ids if not given and if tax ids given in the first
    place or found search for lineage.

    Incorporates functions "taxonomy_metadata", "taxid_metadata" &
    "return_ncbi_taxid".
    """
    taxonomy_csv = taxonomy_metadata(csv_dataframe)
    taxids_csv = taxid_metadata(csv_dataframe)

    print("\nSearching NCBI for taxonomy...")
    combined_lineage = {}
    lineage_custom = {}
    taxids = {}

    for entry, given_taxid in taxids_csv.items():

        if given_taxid != "":
            # If taxid is provided in metadata csv, then add empty entry to
            # lineage_custom dict and taxid to taxid dict
            lineage_custom[entry] = ""
            taxids[entry] = given_taxid

        else:
            taxonomy = taxonomy_csv[entry]
            if not taxonomy:
                # If no tax info is provided in the csv then use provided
                # searchterm
                if isinstance(searchterm, str):
                    tax_id = return_ncbi_taxid(entry, searchterm, email_address)
                    tax_levels = [searchterm]
                else:
                    print(f"No taxonomy information is provided for entry "
                          f"'{entry}' in the csv-file.")
                    term_to_search = input("What term should be searched in "
                                           "NCBI taxonomy?\n")
                    tax_levels = [term_to_search]
                    tax_id = return_ncbi_taxid(entry, term_to_search,
                                               email_address)
            else:
                tax_levels = taxonomy
                tax_id = ""

            c_lineage = []
            n = 0
            while tax_id == "" and n < len(tax_levels):
                # Loop through tax_levels fetching tax_id from NCBI for each.
                # If no hits are returned for a searchterm, give user option to
                # record it as custom lineage information.
                tax_name = tax_levels[n]
                tax_id = return_ncbi_taxid(entry, tax_name, email_address)
                n += 1
                c_lineage.append(tax_name)

            taxids[entry] = tax_id

            if tax_id == "":
                print(f" - For entry '{entry}' no tax id was found.")
                lineage_custom[entry] = "; ".join(tax_levels)
            else:
                lineage_custom[entry] = "; ".join(c_lineage[0:-1])

        combined_lineage[entry] = [taxids[entry], lineage_custom[entry]]

    return combined_lineage

def rejecting_entries(ncbi_lineage, genbank_dict, csv_df, rejection):
    """Print rejected entries to CSV and GenBank files.
    Return df and gb_dict with accepted entries only.
    """
    rejected = []
    new_entries = []
    csv_df.set_index("contigname", inplace=True)

    # Drop rejected rows from DataFrame
    for record, ncbi_info in ncbi_lineage.items():
        if ncbi_info[1] != "" and rejection == "True":
            rejected.append(record)
            entry = (csv_df.loc[record]).values.tolist()
            entry.insert(0, record)
            new_entries.append(entry)
            csv_df.drop([record], inplace=True)

    if len(rejected):
        # Create new DataFrame of rejected entries and drop them from returned
        # DataFrame
        print("\nPrinting rejected entries to CSV file...")
        cols = ['contigname', 'db_id', 'institution_code', 'collection_code',
                'specimen_id', 'morphospecies', 'species', 'genus', 'subfamily',
                'order', 'traptype', 'dev_stage', 'site', 'locality',
                'subregion', 'country', 'latitude', 'longitude', 'size',
                'feeding_behaviour', 'habitat', 'habitat_stratum', 'authors',
                'library', 'datasubmitter', 'projectname', 'genbank_accession',
                'notes']
        new_dataframe = pd.DataFrame(new_entries, columns=cols)
        del new_dataframe['db_id']
        new_dataframe.to_csv('rejected_metadata.csv', index=False)

        for x in rejected:
            print(f" - Entry '{str(x)}' added to CSV file.")

        # Create GenBank file of rejected entries and drop them from returned
        # genbank_dict
        print("\nPrinting rejected entries to GenBank file...")
        rejected_gb_list = [genbank_dict[x] for x in rejected]
        rejected_gb = open("rejected_entries.gb", 'w')
        SeqIO.write(rejected_gb_list, rejected_gb, "genbank")
        rejected_gb.close()

        for x in rejected:
            print(f" - Entry '{str(x)}' added to GenBank file.")
            del genbank_dict[x]

        print("\nTo keep entries with custom lineage information, run rejected "
              "entries again with flag '-r False'.")
    else:
        print("No entries rejected. To reject entries with custom lineage "
              "information, re-run script with flag '-r True'.")

    return genbank_dict, csv_df

def insert_taxid(ncbi_lineage, genbank_dict):
    """Insert tax id into gb data (returned from "ncbi_taxid").
    """
    gb_taxonomy = taxonomy_from_gb(genbank_dict)

    for record, tax_id in gb_taxonomy.items():
        # If there is a taxid in gb, replace it, otherwise create new field for
        # taxid from ncbi or delete if no tax_id given.
        genbank_record = genbank_dict[record]
        ncbi_info = ncbi_lineage[record]
        ncbi_id = ncbi_info[0]

        if tax_id == "":
            # If no tax_id is given in gb file but is in ncbi_lineage,
            # NCBI taxid will be inserted.
            if ncbi_id != "":
                field_given = 0
                for (index, feature) in enumerate(genbank_record.features):
                    if feature.type == "source":
                        feature.qualifiers['db_xref'] = ['taxon:' + ncbi_id]
                        field_given = 1

                if field_given == 0:
                    # if there is no "source" field, add one and add db_xref to
                    # qualifiers.
                    len_record = len(genbank_record.seq)
                    feature_location = FeatureLocation(0, len_record)
                    new_feature = SeqFeature(feature_location, type='source')
                    new_feature.qualifiers['db_xref'] = ['taxon:' + ncbi_id]
                    genbank_record.features.append(new_feature)

        else:
            # If tax_id is given, insert it into gb record qualifiers
            for (index, feature) in enumerate(genbank_record.features):
                if feature.type == "source":
                    if ncbi_id != "":
                        feature.qualifiers["db_xref"] = ['taxon:' + ncbi_id]
                    else:
                        del feature.qualifiers["db_xref"]

    return genbank_dict

def loadnamevariants(source=None):
    """Generates dict of name variants for each gene.
    """
    variants = {}
    types = {}
    products = {}

    # Identify source
    if source is None:
        url = 'https://raw.githubusercontent.com/tjcreedy/biotools/master/' \
              'gene_name_variants.txt'
        source = urllib.request.urlopen(url)
    else:
        source = open(source, 'r')

    # Read source
    for line in source:
        line = line.decode('utf-8').strip()
        meta, listvars = line.split(':')
        name, annotype, product = meta.split(';')
        listvars = [name, product.upper()] + listvars.split(',')
        types[name] = annotype
        products[name] = product
        for v in listvars:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' ' + annotype.upper()]:
                    variants[v + s] = name

    source.close()

    return variants, types, products

def alter_features(genbank_dict):
    """Edit the features in the genbank entries.
    """
    unidentifiable_features = set()
    variants, types, products = loadnamevariants()

    for gb_record, record in genbank_dict.items():
        for (index, feature) in enumerate(record.features):
            if feature.type.upper() == "CDS":
                keys = feature.qualifiers.keys()
                del_features = []
                for key in keys:
                    # Delete unwanted keys
                    if key not in ["gene", "location", "codon_start",
                                   "transl_table", "label",
                                   "product"]:
                        del_features.append(key)

                for f in del_features:
                    del feature.qualifiers[f]

                nametags = ['gene', 'product', 'label', 'standard_name']

                if any(t in feature.qualifiers.keys() for t in nametags):
                    name = 0
                    for t in nametags:
                        if t in feature.qualifiers.keys():
                            name = feature.qualifiers[t][0].upper()
                            break

                    if name in variants.keys():
                        # Standardise names
                        new_name = variants[name]

                        feature.qualifiers["gene"] = [new_name]
                        feature.qualifiers["label"] = [f"{new_name} "
                                                       f"{types[new_name]}"]
                        feature.qualifiers["product"] = [products[new_name]]

                    else:
                        sys.exit(f"ERROR: Unknown gene name for "
                                 f"'{str(gb_record)}' in CDS features: "
                                 f"'{str(name)}'")
                else:
                    unidentifiable_features.add((feature.type,
                                                 feature.location.start,
                                                 feature.location.end))

        if len(unidentifiable_features):
            sys.stderr.write("\nWARNING: The following sequence entries had "
                             "unidentifiable annotations:\n")
            for unidfeats in unidentifiable_features:
                sys.stderr.write(gb_record + ": " + ', '.join(
                    [f + " " + str(s) + "-" + str(e) for f, s, e in
                     unidfeats]) + "\n")

    return genbank_dict

def add_lineage_df(csv_dataframe, combined_lineage):
    """Add columns with tax_id, custom_ and ncbi_lineage to metadata dataframe.
    """
    df_add = pd.DataFrame.from_dict(combined_lineage, orient='index')
    df_add.columns = ["ncbi_taxon_id", "custom_lineage"]
    csv_dataframe.drop(['species', 'genus', 'subfamily', 'family', 'order',
                        'taxid'], axis=1, inplace=True)
    df = pd.merge(df_add, csv_dataframe, left_index=True, right_index=True)
    df.reset_index(level=0, inplace=True)
    df.rename(columns={"index": "contigname"}, inplace=True)

    return df

def reformat_df_cols(df):
    """Reorder DataFrame columns for ingestion into MySQL database.
    """
    df = df[['contigname', 'db_id', 'institution_code', 'collection_code',
             'specimen_id', 'morphospecies', 'ncbi_taxon_id', 'custom_lineage',
             'traptype', 'dev_stage', 'site', 'locality', 'subregion',
             'country', 'latitude', 'longitude', 'size', 'feeding_behaviour',
             'habitat', 'habitat_stratum', 'authors', 'library',
             'datasubmitter', 'projectname', 'genbank_accession', 'notes',
             'version']]

    return df

def load_ids_to_master(new_ids):
    """Loads new db_ids as primary keys into the master table
    """
    con = mdb.connect(host="localhost", user=db_user, passwd=db_passwd,
                      db=db_name)
    with con:
        cur = con.cursor()
        db_ids = "'), ('".join(new_ids)
        sql = f"INSERT INTO master (db_id) VALUES ('{db_ids}');"
        cur.execute(sql)

    return

def load_gb_dict_into_db(genbank_data):
    """Load genbank_data as a dictionary into the mysql database.
    """
    print("\nLoading genbank entries into the database...")

    server = BioSeqDatabase.open_database(driver=db_driver, user=db_user,
                                          passwd=db_passwd, host=db_host,
                                          db=db_name)
    db = server[namespace]
    count = db.load(genbank_data.values())
    server.commit()

    print(" - %i sequences loaded." % count)

    return

def load_df_into_db(csv_dataframe):
    """Loading pandas dataframe with metadata into the database.
    """
    print("\nLoading metadata into the database...")

    engine = create_engine(mysql_engine, echo=False)
    csv_dataframe.to_sql(name='metadata', if_exists='append', index=False,
                         con=engine)

    print(" - %i entries loaded." % len(csv_dataframe.index))

    return

def sql_cols(table, cols, spec):
    """Determine required tables and columns for MySQL query. Construct columns
    string.
    """
    # Reformat inputs
    if spec is None:
        spec = []
    if cols is None:
        cols = []

    # Compile columns requested in query (from cols and spec)
    all_cols = []
    for s in spec:
        req_data = re.split('!=|>=|<=| IN |=|>|<', s)[0]
        if req_data.startswith('(') and req_data.endswith(')'):
            # If a tuple consiting of multiple columns is requested, split it
            # into its constituent columns and add them to the list.
            # E.g. '(bioentry_id, metadata_id) IN ((1, 1), (2, 2))'
            # --> all_cols = [..., bioentry_id, metadata_id]
            split = []
            for col in req_data[1:-1].split(','):
                if col.count('.') >= 1:
                    split.append(col.split('.')[-1].strip())
                else:
                    split.append(col.strip())
            all_cols.extend(split)
        else:
            all_cols.append(req_data)
    all_cols = list(set(all_cols + cols))

    # Unique cols of each table (shared cols assigned to a prioritised table)
    metadata_cols = ['metadata_id', 'contigname', 'db_id', 'institution_code',
                     'collection_code', 'specimen_id', 'morphospecies',
                     'ncbi_taxon_id', 'custom_lineage', 'traptype', 'dev_stage',
                     'site', 'locality', 'subregion', 'country', 'latitude',
                     'longitude', 'size', 'feeding_behaviour', 'habitat',
                     'habitat_stratum', 'metadata.authors', 'library',
                     'datasubmitter', 'projectname', 'genbank_accession',
                     'notes', 'metadata.version']
    bioentry_cols = ['bioentry_id', 'bioentry.biodatabase_id', 'taxon_id',
                     'name', 'accession', 'identifier', 'division',
                     'description', 'version']
    biosequence_cols = ['biosequence.version', 'length', 'alphabet', 'seq']
    comment_cols = ['comment_id', 'comment.bioentry_id', 'comment_text',
                    'comment.rank']
    taxon_cols = ['taxon.taxon_id', 'taxon.ncbi_taxon_id', 'parent_taxon_id',
                  'node_rank', 'genetic_code', 'mito_genetic_code',
                  'left_value', 'right_value']
    taxon_name_cols = ['taxon_name.taxon_id', 'taxon_name.name', 'name_class']

    # Taxonomic queries
    taxonomy = ['subspecies', 'species', 'genus', 'tribe', 'family', 'order',
                'class', 'phylum', 'kingdom', 'superkingdom', 'taxon']

    # Construct columns dictionary (adding prefixes for table joins)
    # E.g. country -> metadata.country
    cols_dict = {}

    for c in all_cols:
        if c == '*':
            continue
        elif c == 'count':
            continue
        elif c in taxonomy:
            mysql_com = f'metadata.TAXON'
        elif c in metadata_cols:
            mysql_com = f'metadata.{c}'
        elif c in bioentry_cols:
            mysql_com = f'bioentry.{c}'
        elif c in biosequence_cols:
            mysql_com = f'biosequence.{c}'
        elif c in comment_cols:
            mysql_com = f'comment.{c}'
        elif c in taxon_cols:
            mysql_com = f'taxon.{c}'
        elif c in taxon_name_cols:
            mysql_com = f'taxon_name.{c}'
        else:
            sys.exit(f"ERROR: Column '{c}' does not exist in the database.")

        if mysql_com.count('.') > 1:
            mysql_com = mysql_com.split('.', 1)[1]

        cols_dict[c] = mysql_com

    # Construct tables list
    tables = list(filter(None, list(
        set([x.split('.')[0] for x in cols_dict.values()] + [table]))))

    # Construct columns string
    if cols == ['*']:
        if len(tables) == 1:
            cols_string = '*'
        else:
            if table is None:
                cols_string = '*'
            else:
                cols_string = f"{table}.*"
    elif cols == ['count']:
        cols_string = "COUNT(*)"
    else:
        cols_string = ', '.join([cols_dict[x] for x in cols])

    return tables, cols_string, cols_dict, spec

def table_join(start, table_list, main_table, shared_col):
    """Consructs string for MySQL table joins

    Arguments:
        - start - starting table/string. E.g. 'bioentry'
        - table_list - list of tables to be joined to starting table/string
                       e.g. ['biosequence']
        - main_table - parent table of all in list. E.g. 'bioentry'
        - shared_col - column shared by all tables in list on which the join
                       is to be performed. E.g. 'bioentry_id'
    """
    n = 0
    table_string = start
    while n < len(table_list):
        new_join = f" JOIN {table_list[n]} ON {main_table}.{shared_col}=" \
                   f"{table_list[n]}.{shared_col}"
        table_string += new_join
        n += 1

    return table_string

def sql_table(tables):
    """Construct table string for MySQL query
    """
    if len(tables) == 1:
        table_string = tables[0]

    elif len(tables) > 1:
        # Lists of tables sharing columns (For duplicates it must be
        # decided which table the column should be assigned to).
        BIOENTRY_ID = ['bioentry', 'bioentry_dbxref',
                       'bioentry_qualifier_value', 'bioentry_reference',
                       'biosequence', 'comment', 'seqfeature']
        TAXON_ID = ['taxon', 'taxon_name']

        # Split required tables into groups according to shared columns
        bios = list(set(tables) & set(BIOENTRY_ID))
        taxons = list(set(tables) & set(TAXON_ID))

        ##JOIN TABLES
        joins = ["metadata"]

        if len(bios) >= 1:
            # Join bios
            if 'bioentry' in bios:
                bios.remove('bioentry')
            start = " JOIN bioentry ON metadata.db_id=bioentry.name"
            bios_join = table_join(start, bios, 'bioentry', 'bioentry_id')
            joins.append(bios_join)

        if len(taxons) >= 1:
            # Join taxons
            if 'taxon' in taxons:
                taxons.remove('taxon')
            start = " JOIN taxon ON metadata.ncbi_taxon_id=taxon.ncbi_taxon_id"
            taxons_join = table_join(start, taxons, 'taxon', 'taxon_id')
            joins.append(taxons_join)

        table_string = ''.join(joins)

    else:
        sys.exit("ERROR: Cannot construct table. Invalid information provided.")

    return table_string

def sql_spec(cols_dict, spec, spec_type):
    """Construct specification string for MySQL query
    """
    # Reformat specifications
    for s in spec:
        split = re.split('!=|>=|<=| IN |=|>|<', s)
        if split[1].isnumeric() or split[1].startswith('('):
            continue
        else:
            ind = spec.index(s)
            op = re.findall('!=|>=|<=| IN |=|>|<', s)[0]
            spec[ind] = f"{split[0]}{op}'{split[1]}'"

    # Condstruct specifications string
    if len(spec) == 0:
        spec_string = ''
    else:
        specs = []
        for x in spec:
            split = re.split('!=|>=|<=| IN |=|>|<', x)
            # If taxonomy query, then add an additional specification
            # containing searchterm.
            if split[0] in ['subspecies', 'species', 'genus', 'tribe', 'family',
                            'order', 'class', 'phylum', 'kingdom',
                            'superkingdom', 'taxon']:
                specs.append(f"metadata.ncbi_taxon_id IN (SELECT DISTINCT "
                             f"include.ncbi_taxon_id FROM taxon INNER JOIN taxon"
                             f" AS include ON (include.left_value BETWEEN "
                             f"taxon.left_value AND taxon.right_value) WHERE "
                             f"taxon.taxon_id IN (SELECT taxon_id FROM "
                             f"taxon_name WHERE name COLLATE LATIN1_GENERAL_CI "
                             f"LIKE '%{split[1][1:-1]}%'))")
                print('Taxonomy searches may take a few minutes...')
            else:
                # Append table names to column names
                pattern = re.compile("|".join(cols_dict.keys()))
                rep_data = []
                for data in split:
                    new_data = pattern.sub(
                        lambda m: cols_dict[re.escape(m.group(0))], data)
                    rep_data.append(new_data)

                # Append modified specification
                specs.append(
                    re.findall('!=|>=|<=| IN |=|>|<', x)[0].join(rep_data))

        # Join specifications
        if spec_type == "output":
            spec_string = f" WHERE ({') AND ('.join(specs)})"

        if spec_type == "update":
            spec_string = ', '.join(specs)

    return spec_string

def construct_sql_output_query(table, cols, spec):
    """Builds MySQL SELECT statement from table, column, and specification
    strings.
    """
    tables, cols_string, cols_dict, spec = sql_cols(table, cols, spec)

    table_string = sql_table(tables)

    spec_string = sql_spec(cols_dict, spec, "output")

    mysql_command = f"SELECT {cols_string} FROM {table_string}{spec_string};"

    return mysql_command

def construct_sql_update_query(table, update, spec):
    """Builds MySQL UPDATE statement from table, column, and specification
    strings.
    """
    tables, _, cols_dict, _ = sql_cols(None, None, spec + update)

    table_string = sql_table(tables)

    spec_string = sql_spec(cols_dict, spec, "output")

    update_string = sql_spec(cols_dict, update, "update")

    mysql_command = f"UPDATE {table_string} SET {update_string}{spec_string};"

    return mysql_command

def construct_sql_delete_query(spec):
    """Builds MySQL DELETE statement from table, column, and specification
    strings.
    """
    tables, _, cols_dict, _ = sql_cols(None, None, spec)

    table_string = sql_table(tables)

    spec_string = sql_spec(cols_dict, spec, "output")

    mysql_command = f"DELETE FROM {table_string}{spec_string};"

    return mysql_command

def fetch_names(mysql_command):
    """Fetch names and corresponding db_id's from database using MySQL query.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)

    with con:
        cur = con.cursor()
        cur.execute(mysql_command)
        records = cur.fetchall()

    names_dict = {row[0]: row[1] for row in set(records)}

    return names_dict

def fetch_recs(names_dict, _all):
    """Fetches a list of SeqRecords from an input dict of record names/db ids
    """
    recs = {}
    server = BioSeqDatabase.open_database(driver=db_driver, user=db_user,
                                          passwd=db_passwd, host=db_host,
                                          db=db_name)
    db = server[namespace]
    for name, db_id in names_dict.items():
        if _all:
            seq_recs = db.get_Seqs_by_acc(db_id)
            for rec in seq_recs:
                recs[rec.id] = rec
        else:
            bio_version, _ = check_current_version(db_id)
            recs[name] = db.get_Seq_by_ver(f'{db_id}.{bio_version}')

    return recs

def execute_query(mysql_query):
    """Connect to db and execute mysql query
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)

    with con:
        cur = con.cursor()
        cur.execute(mysql_query)

    return

def extract_genes(recs, genes):
    """Extracts genes from SeqRecord objects and writes to dict:
    {gene_name :  list_of_sliced_seqrecords, ...}
    """
    if '*' in genes:
        genes = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2',
                 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

    print(f"Extracting {len(genes)} genes: {', '.join(genes)}...")

    subrecs = {}

    for gene in genes:
        extracted_genes = []
        for idd, record in recs.items():
            for feature in record.features:
                if feature.type.upper() == "CDS" \
                        and 'gene' in feature.qualifiers \
                        and feature.qualifiers['gene'][0].upper() == gene.upper():
                    subrec = feature.location.extract(record)
                    subrec.description = re.sub(', [a-z]+ genome$', '',
                                                record.description)
                    subrec.id = record.id
                    subrec.name = record.name
                    extracted_genes.append(subrec)

        subrecs[gene] = extracted_genes

    return subrecs

def csv_from_sql(mysql_command, csv_name):
    """Extract csv file from SQL database.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    cur = con.cursor()
    cur.execute(mysql_command)

    # Write data to CSV file
    with open(f"{csv_name}.csv", "w", newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        headers = [i[0] for i in cur.description]
        rows = [list(row) for row in cur]

        # Delete surrogate keys
        if 'version' in headers or 'metadata_id' in headers:
            indicies = []
            if 'version' in headers:
                version_index = headers.index('version')
                indicies.append(version_index)
            if 'metadata_id' in headers:
                metadata_id_index = headers.index('metadata_id')
                indicies.append(metadata_id_index)
            for index in indicies:
                del headers[index]
                for row in rows:
                    del row[index]

        csv_writer.writerow(headers)
        csv_writer.writerows(rows)

    cur.close()

    return

def seqfile_from_sql(recs_dict, file_name, frmat):
    """Writes list of SeqRecords to a file of chosen format.
    """
    # Specific genes
    if any(isinstance(x, list) for x in recs_dict.values()):
        for gene in recs_dict.keys():
            SeqIO.write(recs_dict[gene], f"{file_name}_{gene}.{frmat}", frmat)

    # Full genome
    else:
        SeqIO.write(recs_dict.values(), f"{file_name}.{frmat}", frmat)

    return

def return_count(mysql_command):
    """Return number of records in the database satisfying a user-supplied
    specification.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    with con:
        cur = con.cursor()
        cur.execute(mysql_command)
        for row in cur:
            print(row[0])

    return

def update_data(metadata, gb_dict):
    """Load updated records into the database.
    """
    ##UPDATE GENETIC DATA
    if gb_dict:

        for rec in gb_dict.values():
            # Update accession and gi
            bio_version, _ = check_latest_version(rec.name)
            bio_version += 1
            rec.id = f"{rec.name}.{bio_version}"
            if 'gi' in rec.annotations.keys():
                del rec.annotations['gi']

        # Load into db
        load_gb_dict_into_db(gb_dict)

    ##UPDATE METADATA
    if metadata is not None:

        metadata.set_index('db_id', inplace=True)

        for db_id in metadata.index:
            # Update version
            _, meta_version = check_latest_version(db_id)
            metadata.ix[db_id, 'version'] = meta_version + 1

        # Load into db
        metadata.reset_index(inplace=True)
        metadata = reformat_df_cols(metadata)
        load_df_into_db(metadata)

    return

def update_master_table(gb_ids, meta_ids, action):
    """Load primary keys (bioentry_id, metadata_id) of most up-to-date versions
    of each record into master table.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    server = BioSeqDatabase.open_database(driver=db_driver, user=db_user,
                                          passwd=db_passwd, host=db_host,
                                          db=db_name)
    db = server[namespace]
    adaptor = db.adaptor

    if action == 'ingest':
        # If the action being performed is data ingestion, then there will only
        # be a single record corresponding to each db_id, so we can just load
        # its primary keys into the master table.
        for db_id in gb_ids:
            bio_id = adaptor.fetch_seqid_by_display_id(1, db_id)
            sql = f"SELECT metadata_id FROM metadata WHERE db_id='{db_id}';"
            with con:
                cur = con.cursor()
                cur.execute(sql)
                meta_id = cur.fetchone()[0]
                sql = f"UPDATE master SET bioentry_id={bio_id}, " \
                      f"metadata_id={meta_id} WHERE db_id='{db_id}';"
                cur.execute(sql)

    else:
        # If the action being performed is an update to records already existing
        # in the database, then load the primary keys of the *latest* version of
        # each into the master table.
        if gb_ids:

            for db_id in gb_ids:
                # Find bioentry_id for latest version
                bio_ver, _ = check_latest_version(db_id)
                bio_id = adaptor.fetch_seqid_by_version(1, f"{db_id}.{bio_ver}")

                # Update master table
                sql = f"UPDATE master SET bioentry_id={bio_id} WHERE " \
                      f"db_id='{db_id}';"
                with con:
                    cur = con.cursor()
                    cur.execute(sql)

        if meta_ids:

            for db_id in meta_ids:
                # Find metadata_id for latest version
                _, meta_ver = check_latest_version(db_id)
                sql = f"SELECT metadata_id FROM metadata WHERE " \
                      f"(db_id='{db_id}') AND (version={meta_ver});"

                with con:
                    cur = con.cursor()
                    cur.execute(sql)
                    meta_id = cur.fetchone()[0]
                    sql = f"""UPDATE master SET metadata_id={meta_id} WHERE 
                        db_id='{db_id}';"""
                    cur.execute(sql)

    return

def rollback_versions(versions_dict):
    """Rollback current version of selected records highlighted in master
    table to a previous version of the users choice.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    server = BioSeqDatabase.open_database(driver=db_driver, user=db_user,
                                          passwd=db_passwd, host=db_host,
                                          db=db_name)
    db = server[namespace]
    adaptor = db.adaptor

    for db_id in versions_dict.keys():
        target_bio_ver = versions_dict[db_id]['b']
        target_meta_ver = versions_dict[db_id]['m']

        # bioentry_id
        if target_bio_ver is not None:
            bio_id = adaptor.fetch_seqid_by_version(1,
                                                    f"{db_id}.{target_bio_ver}")
            update_master = f"""UPDATE master SET bioentry_id={bio_id} 
                WHERE db_id='{db_id}';"""
            with con:
                cur = con.cursor()
                cur.execute(update_master)

        # metadata_id
        if target_meta_ver is not None:
            fetch_id = f"""SELECT metadata_id FROM metadata WHERE (db_id='{db_id}') 
                AND (version={target_meta_ver});"""
            with con:
                cur = con.cursor()
                cur.execute(fetch_id)
                meta_id = cur.fetchone()[0]
                update_master = f"""UPDATE master SET metadata_id={meta_id} 
                    WHERE db_id='{db_id}';"""
                cur.execute(update_master)

    return

def remove_recs(db_ids):
    """Removes all data corresponding to each id in provided ids list from
    database.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    sql = f"DELETE FROM master WHERE db_id IN {tuple(db_ids)};"
    # Deletion from master table automatically cascades to all other tables, as
    # per constraints written into SQL database schema.

    with con:
        cur = con.cursor()
        cur.execute(sql)

    return

def remove_versions(versions_dict):
    """Removes selected versions of records corresponding to user-supplied ids
    from database.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    server = BioSeqDatabase.open_database(driver=db_driver, user=db_user,
                                          passwd=db_passwd, host=db_host,
                                          db=db_name)
    db = server[namespace]
    adaptor = db.adaptor

    for db_id in versions_dict.keys():
        current_bio_ver, current_meta_ver = check_current_version(db_id)
        target_bio_ver = versions_dict[db_id]['b']
        target_meta_ver = versions_dict[db_id]['m']

        if target_bio_ver:
            # If user has specified a bioentry record they wish to delete...
            del_bio_id = adaptor.fetch_seqid_by_version(1,
                                                        f"{db_id}.{target_bio_ver}")

            if target_bio_ver == current_bio_ver:
                # If the target bioentry record is the current for that
                # particular db_id, notify the user that deletion will
                # automatically rollback current version to the *latest*
                # ingested.
                x = input(f"Version {target_bio_ver} is the current version of "
                          f"bioentry '{db_id}'. If removed, the current "
                          f"version will be assumed to be the latest in the "
                          f"database. Do you wish to continue? 'Y'/'N'").upper()
                while not (x == 'Y' or x == 'N'):
                    x = input(f"Would you like to delete the current version "
                              f"of bioentry '{db_id}' ('Y') or cancel the "
                              f"operation ('N')?").upper()
                if x == 'N':
                    sys.exit('Operation cancelled.')
                else:
                    sql_del = f"DELETE FROM bioentry WHERE " \
                              f"bioentry_id={del_bio_id};"
                    with con:
                        cur = con.cursor()
                        cur.execute(sql_del)
                        update_master_table([db_id], None, 'remove')
            else:
                sql_del = f"DELETE FROM bioentry WHERE bioentry_id={del_bio_id};"
                with con:
                    cur = con.cursor()
                    cur.execute(sql_del)

        if target_meta_ver:
            # If user has specified a metadata record they wish to delete...
            sql = f"SELECT metadata_id FROM metadata WHERE (db_id='{db_id}') " \
                  f"AND (version={target_meta_ver});"
            with con:
                cur = con.cursor()
                cur.execute(sql)
                del_meta_id = cur.fetchone()[0]

            if target_meta_ver == current_meta_ver:
                # If the target metadata record is the current for that
                # particular db_id, notify the user that deletion will
                # automatically rollback current version to the *latest*
                # ingested.
                x = input(f"Version {target_meta_ver} is the current version "
                          f"of meta-entry '{db_id}'. If removed, the current "
                          f"version will be assumed to be the latest in the "
                          f"database. Do you wish to continue? 'Y'/'N'").upper()
                while not (x == 'Y' or x == 'N'):
                    x = input(f"Would you like to delete the current version "
                              f"of meta-entry '{db_id}' ('Y') or cancel the "
                              f"operation ('N')?").upper()
                if x == 'N':
                    sys.exit('Operation cancelled.')
                else:
                    sql_del = f"DELETE FROM metadata WHERE " \
                              f"metadata_id={del_meta_id};"
                    with con:
                        cur = con.cursor()
                        cur.execute(sql_del)
                        update_master_table(None, [db_id], 'remove')
            else:
                sql_del = f"DELETE FROM metadata WHERE " \
                          f"metadata_id={del_meta_id};"
                with con:
                    cur = con.cursor()
                    cur.execute(sql_del)

    return