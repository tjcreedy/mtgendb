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

def check_login_details(db_un, dp_pw):
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
        expected_header = ['name', 'specimen', 'morphospecies', 'species',
                           'subfamily', 'family', 'order', 'taxid',
                           'collectionmethod', 'lifestage', 'site', 'locality',
                           'subregion', 'country', 'latitude', 'longitude',
                           'size', 'habitat', 'feeding_behaviour', 'locomotion',
                           'authors', 'library', 'datasubmitter', 'projectname',
                           'accession', 'uin', 'notes']
    else:
        expected_header = ['name', 'db_id', 'morphospecies', 'taxon_id',
                           'custom_lineage', 'specimen', 'collectionmethod',
                           'lifestage', 'site', 'locality', 'subregion',
                           'country', 'latitude', 'longitude', 'size',
                           'habitat', 'feeding_behaviour', 'locomotion',
                           'authors', 'library', 'datasubmitter', 'projectname',
                           'accession', 'uin', 'notes']

    if expected_header != csv_header:
        print("Incorrect header in CSV file.\n")
        sys.exit("Current header is: " + str(csv_header) +
                 "\n\nIt must be as follows: " + str(expected_header))

    return

def matching_inputids(csv_df, gb_dict, action):
    """Check if the GenBank and CSV metadata files have matching entries, and
    that there are no duplicates in the CSV file — duplicates in the GenBank
    file would have already been flagged in gb_into_dictionary().
    """
    new_gb_dict, new_csv_df = [gb_dict, csv_df]

    if action == 'ingest':
        ids_csv = csv_df['name'].values.tolist()
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
            new_csv_df.set_index('name', inplace=True)
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
        mysql_command = f"SELECT * FROM metadata WHERE (name='{idd}') OR " \
                        f"(db_id='{idd}');"
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
        sql_1 = f"SELECT * FROM metadata WHERE name='{acc}';"
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
            # {NC_001: 0, MX432: 0, ... }
        grpn += 1

    groupdups = defaultdict(list)
    for name, group in dupgroups.items():
        groupdups[group].append(name)
        # dups = {1: ['NC_001', 'MX432'], 2: ['NC_008', 'SDJHDFDSJ', 'DHS543']}

    use, reject = [[], []]
    for names in groupdups.values():
        if len(names) == 1:
            use.append(names[0])
        else:
            # If any begin with 'NC_' then use the first of those, otherwise
            # use any
            nc_i = [n for n in names if n.startswith('NC_')]
            if len(nc_i) > 0:
                use.append(nc_i[0])
                rej = [n for n in names if n != nc_i[0]]
            else:
                use.append(names[0])
                rej = names[1:]

            reject.extend(rej)

    if reject:
        print(f"WARNING: The following duplicates have been dropped from your"
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
        print(f"WARNING: The following sequences already exist in the database:"
              f"\n{', '.join(reject)}")

    new_recs = {name: rec for name, rec in records.items() if name not in reject}

    return new_recs, reject

def load_reject_table(rejected):
    """Loads a list of rejected accession numbers into the rejected_accs SQL
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
     - sql - SQL query: "SELECT metadata.name, metadata.db_id FROM... ;"
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
    dict_df = pd.DataFrame(list(dict_new_ids.items()), columns=["name", "db_id"])
    new_csv_df = pd.merge(dict_df, csv_df, on='name')

    return new_csv_df

def taxonomy_metadata(csv_df):
    """Function returning a dictionary with all the taxonomic information for
    each id from metadata ([] if no info given).
    """
    csv_taxa = {}
    ids_csv = csv_df['name'].values.tolist()
    df = csv_df.set_index("name")

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
    gb_met_df.columns = ["name", "taxon_id", "country", "subregion",
                         "collectionmethod", "accession"]

    # Add other metadata columns required by database and fill them with blanks
    # (None/NaN objects)
    for label in ['morphospecies', 'custom_lineage', 'specimen', 'lifestage',
                  'site', 'locality', 'authors', 'size', 'habitat',
                  'feeding_behaviour', 'locomotion', 'library', 'datasubmitter',
                  'projectname', 'uin', 'notes']:
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
    ids_csv = csv_dataframe['name'].values.tolist()
    df = csv_dataframe.set_index("name")

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
    record = handle.read()  # or equivalently: 'record = handle.read()'?
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
                  f"'{entry}'.\n...")

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
    record = handle.read()
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
    csv_df.set_index("name", inplace=True)

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
        cols = ['name', 'db_id', 'specimen', 'morphospecies', 'species',
                'subfamily', 'family', 'order', 'taxid', 'collectionmethod',
                'lifestage', 'site', 'locality', 'subregion', 'country',
                'latitude', 'longitude', 'size', 'habitat', 'feeding_behaviour',
                'locomotion', 'authors', 'library', 'datasubmitter',
                'projectname', 'accession', 'uin', 'notes']
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
                        # TODO: WHY JUST THESE^? WHY NOT 'translation', 'protein_id' etc.?
                        del_features.append(key)

                for f in del_features:
                    del feature.qualifiers[f]

                nametags = ['gene', 'product', 'label', 'standard_name']
                # TODO: HOW COULD THERE POSSIBLY BE 'standard_name' WHEN IT WAS DROPPED IN THE LAST STEP?

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
    df_add.columns = ["taxon_id", "custom_lineage"]
    csv_dataframe.drop(['species', 'subfamily', 'family', 'order', 'taxid'],
                       axis=1, inplace=True)
    df = pd.merge(df_add, csv_dataframe, left_index=True, right_index=True)
    df.reset_index(level=0, inplace=True)
    df.rename(columns={"index": "name"}, inplace=True)

    return df

def reformat_df_cols(df):
    """Reorder DataFrame columns for ingestion into MySQL database.
    """
    df = df[['name', 'db_id', 'morphospecies', 'taxon_id', 'custom_lineage',
             'specimen', 'collectionmethod', 'lifestage', 'site', 'locality',
             'subregion', 'country', 'latitude', 'longitude', 'authors', 'size',
             'habitat', 'feeding_behaviour', 'locomotion', 'library',
             'datasubmitter', 'projectname', 'accession', 'uin', 'notes',
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
        req_data = re.split('=|!=| IN |>|<', s)[0]
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
    # WHICH OF THESE ARE RELEVANT?
    metadata_cols = ['metadata_id', 'name', 'db_id', 'morphospecies',
                     'taxon_id', 'custom_lineage', 'specimen',
                     'collectionmethod', 'lifestage', 'site', 'locality',
                     'subregion', 'country', 'latitude', 'longitude',
                     'metadata.authors', 'size', 'habitat', 'feeding_behavior',
                     'locomotion', 'library', 'datasubmitter', 'projectname',
                     'accession', 'uin', 'notes', 'metadata.version']
    bioentry_cols = ['bioentry_id', 'bioentry.biodatabase_id',
                     'bioentry.taxon_id', 'bioentry.name', 'bioentry.accession',
                     'identifier', 'division', 'description', 'version']
    bioentry_dbxref_cols = ['bioentry_dbxref.bioentry_id',
                            'bioentry_dbxref.dbxref_id', 'bioentry_dbxref.rank']
    bioentry_qualifier_value_cols = ['bioentry_qualifier_value.bioentry_id',
                                     'bioentry_qualifier_value.term_id',
                                     'value', 'bioentry_qualifier_value.rank']
    bioentry_reference_cols = ['bioentry_reference.bioentry_id',
                               'bioentry_reference.reference_id',
                               'bioentry_reference.start_pos',
                               'bioentry_reference.end_pos',
                               'bioentry_reference.rank']
    biosequence_cols = ['biosequence.version', 'length', 'alphabet', 'seq']
    comment_cols = ['comment_id', 'comment.bioentry_id', 'comment_text',
                    'comment.rank']
    seqfeature_cols = ['seqfeature_id', 'seqfeature.bioentry_id',
                       'type_term_id', 'source_term_id', 'display_name',
                       'seqfeature.rank']
    taxon_cols = ['taxon.taxon_id', 'ncbi_taxon_id', 'parent_taxon_id',
                  'node_rank', 'genetic_code', 'mito_genetic_code',
                  'left_value', 'right_value']
    taxon_name_cols = ['taxon_name.taxon_id', 'taxon_name.name', 'name_class']

    # Special queries
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
        elif c in bioentry_cols:  # -name  -taxon_id
            mysql_com = f'bioentry.{c}'
        elif c in bioentry_dbxref_cols:
            mysql_com = f'bioentry_dbxref.{c}'
        elif c in bioentry_qualifier_value_cols:
            mysql_com = f'bioentry_qualifier_value.{c}'
        elif c in bioentry_reference_cols:
            mysql_com = f'bioentry_reference.{c}'
        elif c in biosequence_cols:
            mysql_com = f'biosequence.{c}'
        elif c in comment_cols:
            mysql_com = f'comment.{c}'
        elif c in seqfeature_cols:
            mysql_com = f'seqfeature.{c}'
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
            start = " JOIN taxon ON metadata.taxon_id=taxon.ncbi_taxon_id"
            taxons_join = table_join(start, taxons, 'taxon', 'taxon_id')
            joins.append(taxons_join)

        table_string = ''.join(joins)

    else:
        sys.exit("ERROR: Cannot construct table. Invalid information provided.")

    return table_string

def sql_spec(tables, cols_dict, spec, spec_type):
    """Construct specification string for MySQL query
    """
    """"
    new_spec = []
    for s in spec:
        split = re.split('=|!=| IN |>|<', s)
        if split[1].isnumeric() or split[1].startswith('('):
    """
    spec = [s if re.split('=|!=| IN |>|<', s)[1].isnumeric() or re.split('=|!=| IN |>|<', s)[1].startswith('(') else f"{re.split('=|!=| IN |>|<', s)[0]}{re.findall('=|!=| IN |>|<', s)[0]}'{re.split('=|!=| IN |>|<', s)[1]}'" for s in spec]

    if len(spec) == 0:
        spec_string = ''
    else:
        specs = []
        for x in spec:
            split = re.split('=|!=| IN |>|<', x)
            # If taxonomy query, add additional specification
            if split[0] in ['subspecies', 'species', 'genus', 'tribe', 'family',
                            'order',
                            'class', 'phylum', 'kingdom', 'superkingdom',
                            'taxon']:
                specs.append(
                    f"""metadata.taxon_id IN (SELECT DISTINCT include.ncbi_taxon_id 
                            FROM taxon INNER JOIN taxon AS include ON 
                            (include.left_value BETWEEN taxon.left_value AND taxon.right_value) 
                            WHERE taxon.taxon_id IN (SELECT taxon_id FROM taxon_name WHERE name 
                            COLLATE LATIN1_GENERAL_CI LIKE '%{split[1][1:-1]}%'))""")
                print('Taxonomy searches may take a few minutes...')
            else:
                # Add table names to column names
                pattern = re.compile("|".join(cols_dict.keys()))
                rep_data = []
                for data in split:
                    new_data = pattern.sub(
                        lambda m: cols_dict[re.escape(m.group(0))], data)
                    rep_data.append(new_data)

                # Append modified specification
                specs.append(re.findall('=|!=| IN |>|<', x)[0].join(rep_data))

        # JOIN SPECS
        if spec_type == "output":
            spec_string = f" WHERE ({') AND ('.join(specs)})"

        if spec_type == "update":
            spec_string = ', '.join(specs)

        """
        #PULL LATEST VERSION
        if not _all:

            # Lists of tables with bioentry_id field
            BIOENTRY_ID_TABLES = ['bioentry', 'bioentry_dbxref', 'bioentry_qualifier_value',
                                  'bioentry_reference', 'biosequence', 'comment', 'seqfeature']

            # Determine if any bios tables are being queried
            bios = list(set(tables) & set(BIOENTRY_ID_TABLES))

            if bios:
                current_bio_version = f"(bioentry.bioentry_id, metadata.metadata_id) 
                            IN (SELECT master.bioentry_id, master.metadata_id FROM master 
                            JOIN metadata ON master.metadata_id=metadata.metadata_id 
                            JOIN bioentry ON master.bioentry_id=bioentry.bioentry_id{spec_string});"
                specs.append(current_bio_version)

            else:

                current_meta_version = f"metadata.metadata_id IN (SELECT master.metadata_id 
                            FROM master join metadata on master.metadata_id=metadata.metadata_id{spec_string};)"
                specs.append(current_meta_version)
        """

    return spec_string


def construct_sql_output_query(table, cols, spec):
    # table, cols, spec

    tables, cols_string, cols_dict, spec = sql_cols(table, cols, spec)

    table_string = sql_table(tables)

    spec_string = sql_spec(tables, cols_dict, spec, "output")

    mysql_command = f"SELECT {cols_string} FROM {table_string}{spec_string};"

    return mysql_command


def construct_sql_update_query(table, update, spec):
    # update = ['locomotion=arboreal', 'size=12mm']
    # spec = ['length>25000', 'country=United Kingdom']

    tables, _, cols_dict, _ = sql_cols(None, None, spec + update)

    table_string = sql_table(tables)

    spec_string = sql_spec(tables, cols_dict, spec, "output")

    update_string = sql_spec(tables, cols_dict, update, "update")

    mysql_command = f"UPDATE {table_string} SET {update_string}{spec_string};"

    return mysql_command


def construct_sql_delete_query(spec):
    # update = ['locomotion=arboreal', 'size=12mm']
    # spec = ['length>25000', 'country=United Kingdom']

    tables, _, cols_dict, _ = sql_cols(None, None, spec)

    table_string = sql_table(tables)

    spec_string = sql_spec(tables, cols_dict, spec, "output")

    update_string = sql_spec(tables, cols_dict, update, "update")

    mysql_command = f"DELETE FROM {table_string}{spec_string};"

    return mysql_command


def fetch_names(mysql_command, db_un, db_pw):
    """Fetch names and corresponding db_id's from database using MySQL command
    """
    con = mdb.connect(host=db_host, user=db_un, passwd=db_pw, db=db_name)

    with con:
        cur = con.cursor()
        cur.execute(mysql_command)
        records = cur.fetchall()

    names_dict = {row[0]: row[1] for row in set(records)}

    return names_dict


def fetch_recs(names_dict, db_un, db_pw, _all):
    """Fetches a list of SeqRecords from an input dict of record names/db ids
    """
    # names_dict = {'MH404113': 'GB001', 'KT876913': 'GB007', 'KF364622': 'GB008', 'KT876903': 'GB014'}
    recs = {}

    server = BioSeqDatabase.open_database(driver=db_driver, user=db_un,
                                          passwd=db_pw, host=db_host,
                                          db=db_name)  # driver = "MySQLdb", user = "root", passwd = "mmgdatabase", host = "localhost", db = "mmg_test"
    db = server[namespace]

    for name, db_id in names_dict.items():

        if _all:
            seq_recs = db.get_Seqs_by_acc(db_id)
            for rec in seq_recs:
                recs[rec.id] = rec
        else:
            bio_version, _ = check_current_version(db_id)
            recs[name] = db.get_Seq_by_ver(f'{db_id}.{bio_version}')

    # server.close()

    return recs


def execute_query(mysql_query, db_un, db_pw):
    """Connect to db and execute mysql query
    """
    con = mdb.connect(host=db_host, user=db_un, passwd=db_pw, db=db_name)

    with con:
        cur = con.cursor()
        cur.execute(mysql_query)

    return


def extract_genes(recs, genes):
    """Extracts genes from SeqRecord objects and writes to dict: {gene_name :  list_of_sliced_seqrecords, ...}
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
                if feature.type.upper() == "CDS" and 'gene' in feature.qualifiers and \
                        feature.qualifiers['gene'][
                            0].upper() == gene.upper():  # .upper() not relevant, as alter_features means all in db would be upper, and choice in argparse means all in genes will be upper
                    subrec = feature.location.extract(record)
                    subrec.description = re.sub(', [a-z]+ genome$', '',
                                                record.description)
                    subrec.id = record.id
                    subrec.name = record.name
                    extracted_genes.append(subrec)
        subrecs[gene] = extracted_genes

    return subrecs


def csv_from_sql(mysql_command, csv_name, db_un, db_pw):
    # cols, tablename, csv_name, specs = [None, "metadata", "metadata_output", "subregion='Sabah', collectionmethod='BEATING'"]

    # cur = ((1, 'BIOD00197', 'TEST001', 'Mord33', '219430', '', '1426729', 'YPT', 'Adult', 'Santa Fe', 'PNA', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE2', None, 'Biodiversity Initiative NHM', None, None, '753'), (7, 'BIOD02033', 'TEST007', 'Chrys151', '27439', '', '1721528', 'MALAISE1', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE10', None, 'Biodiversity Initiative NHM', None, None, '766'), (10, 'BIOD01983', 'TEST010', 'Lyc14', '71195', '', '1426866', 'MALAISE1', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE9', None, 'Biodiversity Initiative NHM', None, None, '754'), (25, 'BIOD01740', 'TEST025', 'Ero9', '196992', '', '1721743', 'MALAISE2', 'Adult', 'Cerro Hoya', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE8', None, 'Biodiversity Initiative NHM', None, None, '768'), (26, 'BIOD00733', 'TEST026', 'Lyc6', '71195', '', '1721173', 'MALAISE1', 'Adult', 'Cerro Hoya', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE4', None, 'Biodiversity Initiative NHM', None, None, '762'), (30, 'BIOD01826', 'TEST030', 'Elat5', '30009', '', '1426416', 'FIT2', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE9', None, 'Biodiversity Initiative NHM', None, None, '749'), (31, 'BIOD00778', 'TEST031', 'Curc88', '7042', '', '1426395', 'FIT1', 'Adult', 'Santa Fe', 'P6', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE4', None, 'Biodiversity Initiative NHM', None, None, '749'), (32, 'BIOD00050', 'TEST032', 'Curc36', '7042', '', '1427174', 'FIT2', 'Adult', 'Cerro Hoya', 'P7', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE1', None, 'Biodiversity Initiative NHM', None, None, '757'), (34, 'BIOD00437', 'TEST034', 'Phen3', '94777', '', '1427229', 'FIT1', 'Adult', 'Cerro Hoya', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'H3', None, 'Biodiversity Initiative NHM', None, None, '758'), (35, 'BIOD01314', 'TEST035', 'Mord34', '219430', '', '1426731', 'YPT', 'Adult', 'Santa Fe', 'PNA', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE7', None, 'Biodiversity Initiative NHM', None, None, '753'), (39, 'BIOD02163', 'TEST039', 'Staph127', '29026', '', '1426980', 'MALAISE1', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE10', None, 'Biodiversity Initiative NHM', None, None, '755'), (41, 'BIOD01004', 'TEST041', 'Curc225', '7042', '', '1427202', 'LL1', 'Adult', 'Santa Fe', 'P6', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE5', None, 'Biodiversity Initiative NHM', None, None, '758'), (45, 'BIOD00369', 'TEST045', 'Scar8', '7055', '', '1426136', 'FIT1', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'H3', None, 'Biodiversity Initiative NHM', None, None, '746'), (48, 'BIOD00662', 'TEST048', 'Scar11', '7055', '', '1426272', 'FIT2', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE4', None, 'Biodiversity Initiative NHM', None, None, '748'), (50, 'BIOD01754', 'TEST050', 'Hist17', '110043', '', '1427082', 'FIT2', 'Adult', 'Cerro Hoya', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE8', None, 'Biodiversity Initiative NHM', None, None, '756'), (58, 'BIOD01798', 'TEST058', 'UNK139', '32644', '', '1721589', 'MALAISE1', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE9', None, 'Biodiversity Initiative NHM', None, None, '767'), (59, 'BIOD00057', 'TEST059', 'Elat18', '30009', '', '1721640', 'SLAM1', 'Adult', 'Cerro Hoya', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE1', None, 'Biodiversity Initiative NHM', None, None, '767'), (62, 'BIOD01588', 'TEST062', 'Curc201', '7042', '', '1426914', 'MALAISE2', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE8', None, 'Biodiversity Initiative NHM', None, None, '755'), (63, 'BIOD00890', 'TEST063', 'Chrys74', '27439', '', '1426725', 'YPT', 'Adult', 'Santa Fe', 'PNA', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE5', None, 'Biodiversity Initiative NHM', None, None, '753'), (70, 'BIOD01106', 'TEST070', 'Mel2', '219420', '', '1426972', 'SLAM2', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE6', None, 'Biodiversity Initiative NHM', None, None, '755'), (71, 'BIOD01529', 'TEST071', 'Ero11', '196992', '', '1721537', 'MALAISE2', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE7', None, 'Biodiversity Initiative NHM', None, None, '766'), (73, 'BIOD00313', 'TEST073', 'Elat3', '30009', '', '1426393', 'FIT1', 'Adult', 'Santa Fe', 'P6', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE2', None, 'Biodiversity Initiative NHM', None, None, '749'), (74, 'BIOD00416', 'TEST074', 'Curc312', '7042', '', '1721593', 'MALAISE1', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'H3', None, 'Biodiversity Initiative NHM', None, None, '767'), (77, 'BIOD00870', 'TEST077', 'Chrys102', '27439', '', '1426796', 'MALAISE1', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE5', None, 'Biodiversity Initiative NHM', None, None, '753'), (82, 'BIOD01056', 'TEST082', 'Curc84', '7042', '', '1426375', 'FIT1', 'Adult', 'Santa Fe', 'P5', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE5', None, 'Biodiversity Initiative NHM', None, None, '749'), (85, 'BIOD00217', 'TEST085', 'Curc121', '55867', '', '1721122', 'LL', 'Adult', 'Cerro Hoya', 'P5', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE2', None, 'Biodiversity Initiative NHM', None, None, '762'), (86, 'BIOD01503', 'TEST086', 'Mord26', '219430', '', '1721701', 'MALAISE2', 'Adult', 'Cerro Hoya', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE7', None, 'Biodiversity Initiative NHM', None, None, '768'), (89, 'BIOD01663', 'TEST089', 'Mord63', '219430', '', '1721314', 'SLAM2', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE8', None, 'Biodiversity Initiative NHM', None, None, '764'), (91, 'BIOD02186', 'TEST091', 'Scaph2', '7041', 'Scaphidiidae', '1427161', 'FIT2', 'Adult', 'Cerro Hoya', 'P7', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE10', None, 'Biodiversity Initiative NHM', None, None, '757'), (98, 'BIOD00013', 'TEST098', 'Staph74', '29026', '', '1426424', 'FIT2', 'Adult', 'Santa Fe', 'P1', None, 'Panama', None, None, None, None, None, None, 'D Yeo', 'HE1', None, 'Biodiversity Initiative NHM', None, None, '749'), (104, 'MG193532', 'GB005', None, '2219472', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (113, 'MG193336', 'GB014', None, '2219391', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (115, 'MG193436', 'GB016', None, '2219477', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (116, 'MG193444', 'GB017', None, '2219413', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (125, 'MG193456', 'GB026', None, '2219404', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (145, 'MG193427', 'GB046', None, '2219380', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (149, 'KX035174', 'GB050', None, '1903775', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040002', None, None), (152, 'MG193450', 'GB053', None, '2219467', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (159, 'KX035193', 'GB060', None, '1903793', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040351', None, None), (161, 'KX035187', 'GB062', None, '1903788', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040265', None, None), (164, 'MG193379', 'GB065', None, '2219400', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (169, 'MG193453', 'GB070', None, '2219459', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (174, 'MG193457', 'GB075', None, '2219299', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (181, 'KX035185', 'GB082', None, '1903786', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040174', None, None), (182, 'MG193433', 'GB083', None, '2219360', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (198, 'MG193374', 'GB099', None, '2219449', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (200, 'MG193458', 'GB101', None, '2219406', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (202, 'MG193447', 'GB103', None, '2219465', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (211, 'KX035170', 'GB112', None, '1903771', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1039965', None, None), (216, 'KX035176', 'GB117', None, '1903777', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040049', None, None), (220, 'KX035179', 'GB121', None, '1903780', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040067', None, None), (225, 'MG193439', 'GB126', None, '2219359', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (228, 'MG193530', 'GB129', None, '2219441', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (229, 'KX035186', 'GB130', None, '1903787', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040235', None, None), (237, 'KX035183', 'GB138', None, '1903784', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040118', None, None), (281, 'MG193443', 'GB182', None, '206507', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (282, 'MG193322', 'GB183', None, '2219643', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (283, 'MG193441', 'GB184', None, '2219409', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (288, 'MG193383', 'GB189', None, '2219434', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (295, 'MG193380', 'GB196', None, '2219279', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (298, 'NC_036284', 'GB199', None, '124033', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1043087', None, None), (306, 'KX035184', 'GB207', None, '1903785', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040146', None, None), (310, 'KX035190', 'GB211', None, '1903791', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040331', None, None), (314, 'MG193348', 'GB215', None, '2219369', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (323, 'MG193529', 'GB224', None, '2219280', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (326, 'MG193528', 'GB227', None, '2219405', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (334, 'MG193437', 'GB235', None, '2219461', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (339, 'MG193327', 'GB240', None, '2219347', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (343, 'MG193460', 'GB244', None, '2219418', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (346, 'MG193452', 'GB247', None, '2219476', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (349, 'KX035163', 'GB250', None, '1903764', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1039837', None, None), (354, 'MG193459', 'GB255', None, '2219419', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (355, 'KX035180', 'GB256', None, '1903781', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040069', None, None), (357, 'KX035199', 'GB258', None, '1903798', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1043133', None, None), (368, 'MG193454', 'GB269', None, '2219402', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (369, 'MG193442', 'GB270', None, '2219455', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (374, 'KX035192', 'GB275', None, '1903792', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040341', None, None), (375, 'MG193455', 'GB276', None, '2219446', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (381, 'MG193461', 'GB282', None, '2219374', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (393, 'MG193421', 'GB294', None, '2219372', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (409, 'KX035167', 'GB310', None, '1903768', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1039896', None, None), (425, 'KX035177', 'GB326', None, '1903778', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040052', None, None), (432, 'MG193352', 'GB333', None, '2219398', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (438, 'MG193451', 'GB339', None, '2219470', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (443, 'MG193434', 'GB344', None, '2219462', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (447, 'MG193378', 'GB348', None, '2219401', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (458, 'MG193337', 'GB359', None, '2219441', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (459, 'KX035194', 'GB360', None, '1903794', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1043001', None, None), (461, 'MG193536', 'GB362', None, '2219384', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (465, 'MG193346', 'GB366', None, '2219378', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (466, 'MG193438', 'GB367', None, '2219417', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (467, 'MG193382', 'GB368', None, '2219291', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (479, 'MG193410', 'GB380', None, '2219396', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (482, 'MG193408', 'GB383', None, '2219444', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (487, 'MG193406', 'GB388', None, '2219368', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (488, 'KX035189', 'GB389', None, '1903790', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1040327', None, None), (489, 'MG193446', 'GB390', None, '2219352', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (496, 'MG193445', 'GB397', None, '2219436', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (502, 'MG193531', 'GB403', None, '2219487', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (505, 'MG193448', 'GB406', None, '2219466', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (514, 'KX035195', 'GB415', None, '1903795', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1043031', None, None), (517, 'MG193440', 'GB418', None, '2219393', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (522, 'MG193329', 'GB423', None, '2219350', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (524, 'KX035171', 'GB425', None, '1903772', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1039990', None, None), (528, 'KX035165', 'GB429', None, '1903766', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1039866', None, None), (550, 'MG193426', 'GB451', None, '2219361', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (551, 'KX035172', 'GB452', None, '1903773', None, None, 'Field capture', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, 'BMNH1039994', None, None), (553, 'MG193345', 'GB454', None, '2219373', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (554, 'MG193435', 'GB455', None, '2219453', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (557, 'MG193429', 'GB458', None, '2219451', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (561, 'MG193431', 'GB462', None, '2219392', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (570, 'MG193328', 'GB471', None, '2219452', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None), (576, 'MG193381', 'GB477', None, '2219375', None, None, 'Field trap', None, None, None, None, 'Panama', None, None, None, None, None, None, None, None, None, None, None, None, None))
    # headers = ['metadata_id', 'name', 'db_id', 'morphospecies', 'taxon_id', 'custom_lineage', 'specimen', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'size', 'habitat', 'feeding_behaviour', 'locomotion', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes']
    # Connect to database and execute command
    con = mdb.connect(host="localhost", user=db_un, passwd=db_pw, db=db_name)
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

        csv_writer.writerow(headers)  # Write headers
        csv_writer.writerows(rows)

    cur.close()

    return


def seqfile_from_sql(recs_dict, file_name, frmat):
    # Writes list of SeqRecords to a file of chosen format

    # Specific genes
    if any(isinstance(x, list) for x in recs_dict.values()):
        for gene in recs_dict.keys():
            SeqIO.write(recs_dict[gene], f"{file_name}_{gene}.{frmat}", frmat)

    # Full genome
    else:
        SeqIO.write(recs_dict.values(), f"{file_name}.{frmat}", frmat)

    return


def return_count(mysql_command, db_un, db_pw):
    # mysql_command = "SELECT COUNT(*) FROM biosequence WHERE (length<25000);"

    con = mdb.connect(host="localhost", user=db_un, passwd=db_pw, db=db_name)
    cur = con.cursor()
    cur.execute(mysql_command)

    for row in cur:
        print(row[0])

    cur.close()

    return


def update_data(metadata, gb_dict):
    """Overwrite records in the database
    """
    # 2 cases:
    # 1) Something pulled from db (under db_id), edited, and reingested.
    # 2) A new version of a record (under name) needs to be pushed in - ?????

    ##UPDATE GENETIC DATA
    if gb_dict:

        # Update accession, and gi
        for rec in gb_dict.values():
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

        # Update version
        for db_id in metadata.index:
            _, meta_version = check_latest_version(db_id)
            metadata.ix[db_id, 'version'] = meta_version + 1

        # Load into db
        metadata.reset_index(inplace=True)
        metadata = reformat_df_cols(metadata)
        load_df_into_db(metadata)

    return


def update_master_table(gb_ids, meta_ids, action):
    # Connect to Database
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    # Connect to BioSQL
    server = BioSeqDatabase.open_database(driver=db_driver, user=db_user,
                                          passwd=db_passwd, host=db_host,
                                          db=db_name)  # driver = "MySQLdb", user = "root", passwd = "mmgdatabase", host = "localhost", db = "mmg_test"
    db = server[namespace]
    adaptor = db.adaptor

    if action == 'ingest':
        # If ingest, then just load the primary keys of the record
        # corresponding to each db_id into the master table

        for db_id in gb_ids:
            # Fetch bioentry_id
            bioentry_id = adaptor.fetch_seqid_by_display_id(1, db_id)

            # Fetch metadata_id and update master table
            sql = f"SELECT metadata_id FROM metadata WHERE db_id='{db_id}';"
            with con:
                cur = con.cursor()
                cur.execute(sql)
                metadata_id = cur.fetchone()[0]
                # sql = f"""INSERT INTO master(db_id, bioentry_id, metadata_id)
                #    VALUES ('{db_id}', {bioentry_id}, {metadata_id});"""
                sql = f"""UPDATE master SET bioentry_id={bioentry_id}, 
                    metadata_id={metadata_id} WHERE db_id='{db_id}';"""
                cur.execute(sql)

    else:
        # If update, then load the primary keys of the latest version of each
        # if rollback, check current
        # record into the master table

        if gb_ids:

            for db_id in gb_ids:
                # Find bioentry_id for latest version
                bio_version, _ = check_latest_version(db_id)
                bioentry_id = adaptor.fetch_seqid_by_version(1,
                                                             f"{db_id}.{bio_version}")

                # Update master table
                sql = f"UPDATE master SET bioentry_id={bioentry_id} WHERE db_id='{db_id}';"
                with con:
                    cur = con.cursor()
                    cur.execute(sql)

        if meta_ids:

            for db_id in meta_ids:
                # Find metadata_id for latest version
                _, meta_version = check_latest_version(db_id)
                sql = f"""SELECT metadata_id FROM metadata WHERE (db_id='{db_id}') 
                    AND (version={meta_version});"""

                with con:
                    cur = con.cursor()
                    cur.execute(sql)
                    metadata_id = cur.fetchone()[0]
                    sql = f"""UPDATE master SET metadata_id={metadata_id} WHERE 
                        db_id='{db_id}';"""
                    cur.execute(sql)

    return


"""
def rollback_versions_dict(versions):
    #Rolls back to earlier bioentry version
    
    # versions = {'ME12': {'metadata_version': 4, 'bioentry_version': 2 }, 'ME13': {'metadata_version': 8, 'bioentry_version': 10 }}
    #Bioentry version
    x = input("Would you like to rollback genetic data for selected records? Y/N").upper()
    while not (x == 'Y' or x == 'N'):
        x = input("Type 'Y' if you wish to rollback genetic data for the specified records, or 'N' if you do not.").upper()
    if x == 'Y':
        min_bio_version = min(
            [versions[db_id]['bioentry_version'] for db_id in versions.keys()])
        rollback_no = int(input(f'How many versions would you like to rollback? Max: {min_bio_version}'))

        while rollback_no not in range(min_bio_version + 1):
            rollback_no = int(input(f'ERROR: you can only rollback a maximum of {min_bio_version} versions.'))

        for db_id in versions.keys():
            versions[db_id]['bioentry_version'] -= rollback_no

        print(versions)

    else:
        pass

    #Metadata version
    y = input("Would you like to rollback metadata for selected records? Y/N").upper()
    while not (x == 'Y' or x == 'N'):
        x = input("Type 'Y' if you wish to rollback metadata for the specified records, or 'N' if you do not.").upper()

    if y == 'Y':
        min_meta_version = min(
            [versions[db_id]['metadata_version'] for db_id in versions.keys()])
        rollback_no = int(input(f'How many versions would you like to rollback? Max: {min_meta_version}'))

        while rollback_no not in range(min_meta_version + 1):
            rollback_no = int(input(f'ERROR: You can only rollback a maximum of {min_meta_version} versions'))

        for db_id in versions.keys():
            versions[db_id]['metadata_version'] -= rollback_no

    else:
        pass

    if x == 'N' and y == 'N':
        sys.exit('No rollback selected')

    return
"""


def rollback_versions(versions_dict):
    """Fetch internal ids for set of record names
    1. grab primary keys corresponding to <db_id>.<version_no>
    2. Load into metadata table
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    server = BioSeqDatabase.open_database(driver=db_driver, user=db_user,
                                          passwd=db_passwd, host=db_host,
                                          db=db_name)  # driver = "MySQLdb", user = "root", passwd = "mmgdatabase", host = "localhost", db = "mmg_test"
    db = server[namespace]
    adaptor = db.adaptor

    for db_id in versions_dict.keys():
        target_bio_ver = versions_dict[db_id]['b']
        target_meta_ver = versions_dict[db_id]['m']

        # bioentry_id
        if target_bio_ver is not None:
            bioentry_id = adaptor.fetch_seqid_by_version(1,
                                                         f"{db_id}.{target_bio_ver}")
            update_master = f"""UPDATE master SET bioentry_id={bioentry_id} 
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
                metadata_id = cur.fetchone()[0]
                update_master = f"""UPDATE master SET metadata_id={metadata_id} 
                    WHERE db_id='{db_id}';"""
                cur.execute(update_master)

    return


def remove_recs(db_ids):
    """
    Removes all data corresponding to each id in provided ids list from
    database.
    """
    con = mdb.connect(host=db_host, user=db_user, passwd=db_passwd, db=db_name)
    sql = f"DELETE FROM master WHERE db_id IN {tuple(db_ids)};"

    with con:
        cur = con.cursor()
        cur.execute(sql)

    return


def remove_versions(versions_dict):
    """
    Removes selected versions of records corresponding to provided ids from
    database.
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

            del_bio_id = adaptor.fetch_seqid_by_version(1, f"{db_id}.{target_bio_ver}")

            if target_bio_ver == current_bio_ver:
                x = input(f"""Version {target_bio_ver} is the current version
                    of bioentry '{db_id}'. If removed, the current version will be 
                    assumed to be the latest in the database. Do you wish to 
                    continue? 'Y'/'N'""").upper()
                while not (x == 'Y' or x == 'N'):
                    x = input(f"""Would you like to delete the current version
                        of bioentry '{db_id}' ('Y') or cancel the operation ('N')?""").upper()
                if x == 'N':
                    sys.exit('Operation cancelled.')
                else:
                    sql_del = f"DELETE FROM bioentry WHERE bioentry_id={del_bio_id};"
                    with con:
                        cur = con.cursor()
                        cur.execute(sql_del)
                        update_master_table([db_id], None, 'remove')
                        """
                        latest_bio_ver, _ = check_latest_version(db_id)
                        latest_bio_id = adaptor.fetch_seqid_by_version(1, f"{db_id}.{latest_bio_ver}")
                        sql_update = f"UPDATE master SET bioentry_id={latest_bio_id} 
                            WHERE db_id='{db_id}';"
                        cur.execute(sql_update)
                        """
            else:
                sql_del = f"DELETE FROM bioentry WHERE bioentry_id={del_bio_id};"
                with con:
                    cur = con.cursor()
                    cur.execute(sql_del)

        if target_meta_ver:

            sql = f"""SELECT metadata_id FROM metadata WHERE (db_id='{db_id}') 
                AND (version={target_meta_ver});"""
            with con:
                cur = con.cursor()
                cur.execute(sql)
                del_meta_id = cur.fetchone()[0]

            if target_meta_ver == current_meta_ver:
                x = input(f"""Version {target_meta_ver} is the current version
                    of meta-entry '{db_id}'. If removed, the current version will be 
                    assumed to be the latest in the database. Do you wish to 
                    continue? 'Y'/'N'""").upper()
                while not (x == 'Y' or x == 'N'):
                    x = input(f"""Would you like to delete the current version
                        of meta-entry '{db_id}' ('Y') or cancel the operation ('N')?""").upper()
                if x == 'N':
                    sys.exit('Operation cancelled.')
                else:
                    sql_del = f"DELETE FROM metadata WHERE metadata_id={del_meta_id};"
                    with con:
                        cur = con.cursor()
                        cur.execute(sql_del)
                        update_master_table(None, [db_id], 'remove')
                        """
                        _, latest_meta_ver = check_latest_version(db_id)
                        fetch_id = f"SELECT metadata_id FROM metadata WHERE (db_id='{db_id}') 
                            AND (version={latest_meta_ver});"
                        cur.execute(fetch_id)
                        latest_meta_id = cur.fetchone()[0]
                        sql_update = f"UPDATE master SET metadata_id={latest_meta_id} 
                            WHERE db_id='{db_id}';"
                        cur.execute(sql_update)
                        """

            else:
                sql_del = f"DELETE FROM metadata WHERE metadata_id={del_meta_id};"
                with con:
                    cur = con.cursor()
                    cur.execute(sql_del)

    return
