#!/urs/bin/env python

"""Functions to modify and upload a genbank-file and metadata in a csv-file into a MySQL database.
"""


###IMPORT PACKAGES


import sys, time, urllib.request, csv, re

##Import modules for handeling the genbank files:
from Bio import SeqIO, Entrez, Alphabet                             # Entrez is a molecular biology database system that provides integrated access to NCBI’s databases such as PubMed, GenBank, GEO, and many others.
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

##Import modules for handeling the csv files:
import pandas as pd
import numpy as np

##Import modules for dealing with the mysql database:
import MySQLdb as mdb
from BioSQL import BioSeqDatabase
from sqlalchemy import create_engine, update


###VARIABLES

db_driver = "MySQLdb"                                                # This is to do with talking to the database, which you can think of as a separate
db_passwd = "mmgdatabase"                                            # external program that's always running and available to receive, store and give
db_host = "localhost"                                                # out data. These are the settings for it
db_user = "root"
db_name = "mmg_test"
mysql_engine = "mysql+mysqldb://root:mmgdatabase@localhost/mmg_test"

namespace = "mmg"


###DEFINE FUNCTIONS


def gb_into_dictionary(gb_filename, key):
    """Take a file in genbank-format and load it into a dictionary, define function for keys (default: by name in LOCUS).
    
    The reason you provide a function for finding the key and not just the key itself is that the object 'record' doesn't exist at the point seqIO.to_dict
    is called. It doesn't know that you're talking about the objects currently being parsed. That's what the function does - the parser applies that 
    function to each item (i.e. each entry in the genbank file), calling each item a 'record'. Then it extracts out the description, so that seqio.to_dict 
    uses that as the key for that item
    
    KEY_FUNCTION - Optional callback function which when given a SeqRecord should return a unique key for the dictionary.
    """

    global gb_dictionary
    if key == "LOCUS":

        def get_seqname(record):
            #Given a SeqRecord, return the sequence name as a string.
            seqname = str(record.name)                                                                                       # in a genbank file, ID is the accession number and name is the locus name. str() is used here as the name may contain only numbers, in which case it may get parsed in as an int or float. The str is just to make sure it's the right format for that edge case.
            return seqname

        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"), key_function = get_seqname)

    elif key == "ACCESSION":
        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"))                                                   # 'id' is used here by default (as it is whenever key_function is omitted)

    elif key == "DEFINITION":

        def get_definition(record):
            #Given a SeqRecord, return the definition as a string.
            definition = str(record.description)                                                            #4 Isn't the description already a string? (Says it is on BioPython)
            return definition

        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"), key_function = get_definition)

    return gb_dictionary


def text_to_list(accessions):
    """Converts text-file of IDs (one per line) into a list, exiting if duplicates are present.
    """
    #Create a comma-delimited list of accession numbers from text file (stripping any blank spaces/empty lines).
    with open(accessions, "r") as acc:
        accs = acc.read()
        striplist = lambda lis:[x.strip() for x in lis]
        acc_list = list(filter(None, striplist(accs.split('\n'))))
        #e.g. acc_list = ['MH431903', 'MH281574', 'MH281574', 'KT876903']

    #Check for duplicates
    duplicates = set([idd for idd in acc_list if acc_list.count(idd) > 1])

    if len(duplicates):
        #If there are any duplicates, drop these proceed with unque accessions only.
        print("WARNING: Duplicates detected in text file: '" + "', '".join(duplicates) + "'. Proceeding with unique accessions only.")
        acc_list = list(set(acc_list))

    return acc_list


def check_acc_format(acc_list):
    """Checks each accession in list satisfies genbank format, dropping any that don't.
    """
    Accs = []

    for acc in acc_list:
        z = re.match("[A-Z]{1,3}_?[0-9]{5,8}$", acc)
        if z:
            Accs.append(acc)
        else:
            print(f"WARNING: Accession '{acc}' dropped due to incorrect GenBank format.")

    return Accs


def correct_header(csv_dataframe):
    """Check metadata file has correct columns #by creating a list of its column headers and checking it against a list of the expected headers.
    """
    # csv_dataframe = csv_df

    csv_header = csv_dataframe.columns.values.tolist()                                     # LUKE: creates list of csv column headers.      ## Pandas index.tolist() Converts dataframes column values to list, like below.  ###6 LUKE: So why use the .values function when .columns returns an index already?
    expected_header = ['name', 'specimen', 'morphospecies', 'species', 'subfamily', 'family', 'order', 'taxid', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'size', 'habitat', 'feeding_behaviour', 'locomotion', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes']

    if expected_header != csv_header:
        print("Incorrect header in CSV file.\n")
        sys.exit("Current header is: " + str(csv_header) + "\n\nIt must be as follows: " + str(expected_header))

    return()


def matching_inputids(csv_dataframe, gb_dictionary, action):
    """Check if the GenBank and CSV metadata files have matching entries, and that there are no duplicates in the CSV file. (Duplicates in the GenBank file would
    have generated an error when first loaded into a dictionary in gb_into_dictionary.)
    """
    # gb_dictionary, csv_dataframe = [gb_dict, csv_df]

    new_gb_dict, new_csv_df = [gb_dictionary, csv_dataframe]

    if action == 'ingest':
        ids_csv = csv_dataframe['name'].values.tolist()
    else:
        ids_csv = csv_dataframe['db_id'].values.tolist()

    unique = list(set(ids_csv + list(gb_dictionary.keys())))                        #all IDs in both files (the union).
    shared = list(set(ids_csv) & set(gb_dictionary.keys()))                         #all IDs shared by by both files (the intersection).
    discrepant_ids = [a for a in unique if a not in shared]                         #all IDs that aren't shared by both files.
    csv_duplicates = set([idd for idd in ids_csv if ids_csv.count(idd) > 1])  #all duplicate IDs in the CSV file.

    if len(csv_duplicates):
        # Check IDs are unique in the CSV metadata file.
        sys.exit("ERROR: There are multiple rows sharing the same name in your CSV file: " + ', '.join(csv_duplicates) + ". IDs must be unique.")

    if len(shared) != len(unique):
        # If the IDS in the CSV and GenBank files are not identical...
        x = input("WARNING: Your CSV and GenBank files contain different entries.\nWould you like to ignore these and proceed with shared entries only ('P') or cancel the operation ('C')?\n?>").capitalize()
        while not (x == 'C' or x == 'P'):
            x = input("Type 'P' to ignore discrepant entries and proceed, or type 'C' to cancel the operation.\n?>").capitalize()
        if x == 'C':
            sys.exit("Operation cancelled.")
        else:
            #Return new gb_dict with discrepant entries deleted
            for gb_record in gb_dictionary:
                if gb_record in discrepant_ids:
                    print(f" - Skipping entry '{gb_record}' as it appears in the GenBank file but not the CSV file.")

            new_gb_dict = {key: gb_dictionary[key] for key in gb_dictionary if key not in discrepant_ids}

            #Return new csv_dataframe with discrepant entries deleted
            for contigname in ids_csv:
                if contigname in discrepant_ids:
                    print(f" - Skipping entry '{contigname}' as it appears in the CSV file but not the GenBank file.")

            df_missing_ids = [b for b in discrepant_ids if b not in gb_dictionary]
            new_csv_df.set_index('name', inplace=True)
            new_csv_df.drop(df_missing_ids, inplace=True)   #drop rows with ids that aren't shared by both files #THIS NEEDS TO DROP IDS THAT ARE IN DF AND NOT IN GB_DICT, AS IT CAN'T DROP ITEMS THAT AREN'T IN IT
            new_csv_df.reset_index(inplace=True)

    return new_csv_df, new_gb_dict


def new_ids(genbank_dict, prefix, startvalue, padding):
    """Check if the new ids are looking fine given the input (prefix, startvalue, padding) by the user.
    """

    dict_new_ids = {}

    if len(prefix) > 4:
        sys.exit("ERROR: The prefix (-p) is too long (it should be no more than 4 characters).")

    if prefix.isalpha() == False:
        sys.exit("ERROR: The prefix (-p) should only consist of letters.")

    if str(startvalue).isdigit() == False:
        sys.exit("ERROR: The number (-n) should only consist of digits.")

    if len(str(len(genbank_dict.keys()) + int(startvalue) - 1)) > int(padding):
        # Checks there's not too many sequences to assign serial numbers to given the padding. E.g. If we've set the padding to 3 but there are 4295 sequences in the genbank, we can't generate enough numbers to assign them all new names.             
        sys.exit("ERROR: A 0-padding by " + str(padding) + " digits is not sufficient for the ingestion of " + str(len(genbank_dict)) + " new entries.")
        #sys.exit("ERROR: A 0-padding by %s digits is not sufficient for the ingestion of %d new entries." % (padding, len(genbank_dict)))

    if len(str(startvalue)) > int(padding):
        sys.exit("ERROR: The starting number " + str(startvalue) + " exceeds the digits for 0-padding (" + str(padding) + " digits).")
        #sys.exit("ERROR: The starting number %s exceeds the digits for 0-padding (%s digits)." % (startvalue, padding))

    for gb_record, record in genbank_dict.items():
        record_name = record.name
        pad = "{0:0" + str(padding) + "d}"
        padded = pad.format(int(startvalue))
        new_name = prefix + str(padded)
        dict_new_ids[record_name] = new_name
        startvalue = int(startvalue) + 1

    return dict_new_ids

"""
The condition in 4th if block has been changed as the previous condition - len(str(len(genbank_dict.keys()))) <= int(padding) - breaks down in 
the case where the startvalue is high enough to push the total no. of sequences over the specified padding limit. E.g. imagine one has already ingested 500 
sequences into the database with padding=3 (say, BIOD001 - BIOD500), and wishes to ingest 510 more. There are 999 serial number options for padding=3, so if he 
sets the padding to 3 and the startvalue to 501, there will now be more sequences in the database (1010) than names to assign to them, but an error won't be 
generated as the new IDs still pass the condition len(str(len(genbank_dict.keys())))=3 <= 3=int(padding). This is a problem resolved by the new condition - 
len(str(len(genbank_dict.keys()) + startvalue - 1)) <= int(padding): In our example: len(str(len(genbank_dict.keys()) + startvalue - 1)) =
len(str(510 + 501 -1)) = len(str(1010)) = 4  <= 3 = int(padding). 4 is not smaller or equal to 3, and so our condition is violated and an error is generated
asking the user to increase the padding.
"""


def check_ids(db_un, db_pw, ids_list, action):
    """Check if the database ids already exist in the database.
    """
    for idd in ids_list:
        mysql_command = f"SELECT * FROM metadata WHERE (name='{idd}') OR (db_id='{idd}');"
        con = mdb.connect(host=db_host, user=db_un, passwd=db_pw, db=db_name)

        with con:
            cur = con.cursor()
            cur.execute(mysql_command)
            result = cur.fetchone()

        if action == 'replace':
            if result is None:
                sys.exit(f"ERROR: The name/id '{idd}' is missing from the database.")

        elif action == 'ingest':
            if result is not None:
                sys.exit(f"ERROR: The name/id '{idd}' is already present in the database:\n{str(result)}")

        else:
            sys.exit(f"ERROR: Unrecognized action: '{action}'")

    return


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


def change_ids_genbank(genbank_dict, dict_new_ids, key):
    """Change the ids in the genbank file to the new database ids (LLLLNNNNN).
    """
    # genbank_dict, dict_new_ids, key = [new_dict, dict_new_ids, args.key]

    for gb_record, record in genbank_dict.items():                                                                    # L: for each key in genbank_dict, let 'record' equal the value for that key.

        #Change the name (record.name) (to key) of each entry to new database id (to value)
        record.name = dict_new_ids[record.name]                                                                     # record_name = record.name = dict_new_ids[record_name] = new_name

        #Change the accession (record.id) to the new database id including version number (.0)
        new_accession = record.name + ".0"
        record.id = new_accession
        record.annotations["accessions"] = [record.name]

        #Change the version
        record.annotations["sequence_version"] = 0       # Genbank files include a version annotation. This sets it to 0, because we're going to use that ourselves to record our own versions of an entry.

        if key == "DESCRIPTION":
            #Erase the description as it contained the input id
            record.description = ""

    return


def change_names_csv(csv_dataframe, dict_new_ids):
    """Create dictionary with new database ids and old input names to df.
    """

    dict_df = pd.DataFrame(list(dict_new_ids.items()), columns = ["name", "db_id"])  # list(dict_new_ids.items()) creates a list of tuples of the dict pairs - (key, value)
    #creates a 2-column dataframe with the "name" column featuring all the the old ids and the "db_id" column featuring the new ids.

    #Merge csv_dataframe and dict_df -> append new database ids as a column
    new_csv_df = pd.merge(dict_df, csv_dataframe, on = 'name')                       # this merges both dataframes on the 'name' (old id) column, which in effect simply adds the "db_id" (new id) column in dict_df onto the csv_dataframe

    return new_csv_df  # returns csv_dataframe with the new ids added on in a single new column called "db_id"

def taxonomy_metadata(csv_dataframe):
    """Function returning a dictionary with all the taxonomic information for each id from metadata ([] if no info given).
    """

    csv_taxa = {}
    ids_csv = csv_dataframe['name'].values.tolist()     # ids_csv = list of entries in the csv_dataframe 'name' column (is the values bit necessary? Can only see this being the case if the column entries are dictionaries - but then there would be parentheses after i.e. .values())
    df = csv_dataframe.set_index("name")                # df = csv_dataframe with its index set to the "name" column

    for id in ids_csv:                 #for each name in the list ids_csv...

        tax_list = []
        species = df.loc[id, 'species']          # species = the df entry for the 'id' row in the 'species' column. (i.e. the species of the id)
        subfamily = df.loc[id, 'subfamily']      # subfamily = the df entry for the 'id' row in the 'subfamily' column (i.e. the subfamily of the id)
        family = df.loc[id, 'family']            # the family of the id
        order = df.loc[id, 'order']              # the order of the id

        if not pd.isnull(species):               # if the species is not missing/empty...
            tax_list.append(species)             # add it to the tax_list.
        if not pd.isnull(subfamily):             # if the subfamily is not missing/empty...
            tax_list.append(subfamily)           # add it to the tax_list.
        if not pd.isnull(family):                # if the family is not missing/empty...
            tax_list.append(family)              # add it to the tax_list.
        if not pd.isnull(order):                 # if the order is not missing/empty...
            tax_list.append(order)               # add it to the tax_list.

        csv_taxa[id] = tax_list                  # add a key 'id' to the dictionary csv_taxa with value 'tax_list' (i.e. a list of all that ids taxonomic info: [species, subfamily, family, order]).

    return(csv_taxa)                             # return the dictionary csv_taxa.


def chunker(seq, size):
    """Splits input list/set into subsets of specified size
    """
    if type(seq) is set:
        seq = list(seq)

    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def return_gb_data(acc_list, email):
    """Takes list of of IDs/accessions and returns dict of corresponding GenBank entries.
    """
    #email = args.users_email

    print("Fetching GenBank records...")

    records = {}
    Entrez.email = email

    for accset in chunker(acc_list, 200):
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
    #Extract metadata to separate dict
    gb_metadata = {}
    for acc, record in records.items():
        for feature in record.features:
            if feature.type == "source":

                record_metadata = []

                if 'db_xref' in feature.qualifiers:
                    taxid = feature.qualifiers['db_xref'][0].split(":")[1]
                else:
                    print(f"WARNING: 'db_xref' absent from '{acc}' qualifiers.")
                    taxid = None

                if 'country' in feature.qualifiers:
                    country = feature.qualifiers['country'][0]
                    if ': ' in country:
                        country = feature.qualifiers['country'][0].split(": ")[0]
                        subregion = feature.qualifiers['country'][0].split(": ")[1]
                    else:
                        subregion = None
                else:
                    country, subregion = [None, None]

                if 'isolation_source' in feature.qualifiers:
                    traptype = feature.qualifiers['isolation_source'][0]
                else:
                    traptype = None

                if 'specimen_voucher' in feature.qualifiers:
                    accession = feature.qualifiers['specimen_voucher'][0].translate({ord(c):None for c in ' -.:'})     #Strips all instances of ' ', '-', '.' and ':' from string.
                else:
                    accession = None

                record_metadata.append(taxid)
                record_metadata.append(country)
                record_metadata.append(subregion)
                record_metadata.append(traptype)
                record_metadata.append(accession)

                gb_metadata[acc] = record_metadata

    #Write to DataFrame
    gb_met_df = pd.DataFrame.from_dict(gb_metadata, orient='index')
    gb_met_df.reset_index(inplace=True)
    gb_met_df.columns = ["name", "taxon_id", "country", "subregion", "collectionmethod", "accession"]

    #Add other metadata columns required by database and fill them with blanks (None/NaN objects)
    for label in ['morphospecies', 'custom_lineage', 'specimen', 'lifestage', 'site', 'locality', 'authors', 'size', 'habitat', 'feeding_behavior', 'locomotion', 'library', 'datasubmitter', 'projectname', 'uin', 'notes']:
        gb_met_df[label] = None

    for header in ['latitude', 'longitude']:
        gb_met_df[header] = np.nan

    return gb_met_df


def taxonomy_from_gb(genbank_dict):
    """Function returning dictionary of taxids from genbank file if it has it (or "" if not).
    """

    gb_taxa = {}

    for gb_record, record in genbank_dict.items():

        tax_id = ""

        for (index, feature) in enumerate(record.features):                                                                                # Need to understand more about record attributes and their functions for this. (Features are SeqFeature objects)

            if feature.type == "source":        # .type - the specified type of the feature (ie. CDS, exon, repeat...)

                if "db_xref" in feature.qualifiers:        # These are analogous to the qualifiers from a GenBank feature table. The keys of the dictionary are qualifier names, the values are the qualifier values.
                    xref = feature.qualifiers["db_xref"]  # if key 'db_xref' is in dictionary 'feature.qualifiers', 'xref' = its value. xref is a list of strings: ["taxon:nnnnnnn", "taxon2:mmmmmmm"...]

                    if len(xref) > 1:                                                            # if length of 'xref' > 1...
                        print("There are several ids in the db_xref 'source' information.")
                    else:
                        tax_id = xref[0].split(":")[1]   # xref is a list of 1 string: ["taxon:nnnnnnn"], so tax_id= "nnnnnn"

        gb_taxa[record.name] = tax_id          # gb_taxa = {record.name: tax_id, ...}

    return gb_taxa


def taxid_metadata(csv_dataframe):
    """Return dictionary with taxids provided in metadata table ("" if not given).
    """

    csv_taxids = {}
    ids_csv = csv_dataframe['name'].values.tolist()  # ids_csv = list of all entries in csv_dataframe column "name".
    df = csv_dataframe.set_index("name")             # df = csv_dataframe with "name" column as index

    for id in ids_csv:                               # for each id in ids_csv...

        taxid = df.loc[id, 'taxid']                  # taxid = entry in 'id' row & 'taxid' column of dataframe 'df'

        if not pd.isnull(taxid):                     # if the taxid entry is not empty...
            csv_taxids[id] = taxid                   # then add 'id':'taxid' to the (currently empty) dictionary csv_taxids.

        else:                                        # if the taxid is empty...
            csv_taxids[id] = ""                      # then add 'id':'' to csv_taxids.

    return csv_taxids


def return_ncbi_taxid(entry, searchterm, email_address):
    """For each entry get tax_id from NCBI taxonomy based on taxonomic information.
    """

    Entrez.email = email_address                                                                     # optional email parameter so the NCBI can contact you if there is a problem. Good practise to always tell NCBI who you are.
    handle = Entrez.esearch(db = "taxonomy", retmax = 2, term = searchterm)     # what's retmax?     # Search in Taxonomy database for publications relating to 'searchterm'. [.esearch() searches the Entrez databases.]
    record = Entrez.read(handle)                                                                     # or equivalently: 'record = handle.read()'?
    handle.close()

    id_list = record["IdList"]            # id_list = list of ids (e.g. ’28011774’) obtained from the search.

    if len(id_list) == 0:                                                                            # if the search found nothing...
        # Give user option to use unsuccessful searchterm as custom lineage info or cancel the operation.
        x = input(f" - No hits found for search term '{searchterm}' in NCBI taxonomy.\n   Would you like to record this as custom lineage information for entry '{entry}' and proceed ('P') or cancel the operation ('C')?\n ?>").capitalize()
        while not (x == 'P' or x == 'C'):
            x = input(f"   Type 'P' to record '{searchterm}' as custom lineage information for entry '{entry}' or 'C' to cancel the operation.\n ?>").capitalize()
        if x == 'C':
            sys.exit("\nOperation cancelled.")
        else:
            tax_id = ""
            print(f" - '{searchterm}' saved to custom lineage information for '{entry}'.\n...")

    elif len(id_list) > 1:                                                                           # elif the search found more than one result...
        print(f" - Multiple hits found for search term '{searchterm}' in NCBI taxonomy.")
        tax_id = ""                                                                                  # set tax_id to nothing

    elif len(id_list) == 1:                                                                          # elif the search found just one result...
        tax_id = id_list[0]                                                                          # then set the tax_id to the one id obtained from the search.

    time.sleep(0.5)                                                                                  # Pause python program for 0.5 secs

    return tax_id                                                                                    # return tax_id obtained from search


def return_ncbi_lineage(searchterm, email_address):
    """Search NCBI for lineage information given a tax id.
    """

    Entrez.email = email_address                                                             # Always tell NCBI who you are.

    handle = Entrez.efetch(db = "taxonomy", id = searchterm)                    # retrieve the full specified record from the database "taxonomy".
    record = Entrez.read(handle)                                                             # or equivalently: 'record = handle.read()'?
    handle.close()

    if len(record) == 0:                                                                     # if the record is empty...
        print(f"No hits found for tax id '{searchterm}' in NCBI taxonomy.")

    elif len(record) > 1:
        print(f"Multiple hits found for tax id '{searchterm}' in NCBI taxonomy.")

    elif len(record) == 1:                                                                   # if the search found 1 hit...
        taxonomy = record[0]["Lineage"]                                                      # let taxonomy = value for key "Lineage" in first element (dict/list?) of the record for that hit. #?? Is record a dictionary/list of nested dictionaries?
        taxon = record[0]["ScientificName"] #Add the taxon itself to the lineage             # let taxon = [value for key "ScientificName" in first element (dict/list?) of record.
        lineage = taxonomy + "; " + taxon                                                    # lineage = "taxonomy; taxon"

    time.sleep(0.5)                                                                          # Pause python program for 0.5 secs

    return lineage


def get_ncbi_lineage(csv_dataframe, email_address, searchterm):
    """Search on NCBI for tax ids if not given and if tax ids given in the first place or found search for lineage.

    Incorporates functions "taxonomy_metadata", "taxid_metadata" & "return_ncbi_taxid".
    """

    taxonomy_csv = taxonomy_metadata(csv_dataframe) #Get all tax info from csv as dict
    taxids_csv = taxid_metadata(csv_dataframe)      #Get all tax ids from csv as dict

    print("\nSearching NCBI for taxonomy...")
    combined_lineage = {}
    #lineage_ncbi = {}
    lineage_custom = {}
    taxids = {}

    for entry, given_taxid in taxids_csv.items():                                       # for each key (id) in the dictionary of taxids...                               # set given_taxid = the value for that key (the taxid info) in taxids_csv

        if given_taxid != "":                                          # if taxid is provided in metadata csv...
            #Search given tax id on ncbi & add ncbi_lineage
            #ncbi_l = return_ncbi_lineage(given_taxid, email_address)       # ncbi_1 = lineage for given_taxid = "taxonomy; taxon"   # [need to understand this function better]
            #lineage_ncbi[entry] = ncbi_l                                   # Add key-value pair to empty lineage_ncbi dict:    (entry; lineage )    i.e. (id: "taxonomy; taxon")
            lineage_custom[entry] = ""                                     # Add key-value pair to empty lineage_custom dict:  (entry; "")          i.e. (id: "")
            taxids[entry] = given_taxid                                    # Add key-value pair to empty taxids dict:          (entry; given_taxid) i.e. (id: given_taxid)   ## HOW IS THIS ANY DIFFERENT FROM TAXIDS_CSV, AT THIS POINT?

        else:                                                          # if no taxid in metadata csv...

            taxonomy = taxonomy_csv[entry]                                 # taxonomy = list of tax levels from taxonomy_csv dict (or an empty list [])

            if taxonomy == []:                                             # if no tax info is provided in the csv... (then reject or use user input!)

                if isinstance(searchterm, str):                                # if the provided searchterm is a string...
                    tax_id = return_ncbi_taxid(entry, searchterm, email_address)      # tax_id = id number returned from 'return_ncbi_taxid' function given searchterm.
                    tax_levels = [searchterm]                                  # Set tax_levels = [specified searchterm]

                else:                                                          # if the provided searchterm is not a string...
                    print(f"For entry '{entry}' no information about the taxonomy is given in the csv-file.")
                    term_to_search = input("What term should be searched in NCBI taxonomy?\n")       # User input
                    tax_levels = [term_to_search]                                                    # Set tax_levels = User input
                    tax_id = return_ncbi_taxid(entry, term_to_search, email_address)                        # Put user input and email through "return_ncbi_taxid" function to obtain tax_id number. ##If nothing is found the tax id is ""

            else:                                                          # if taxonomy info *is* provided in list in the taxonomy_csv dict...
                #Tax info is searched on ncbi
                tax_levels = taxonomy                                      # Set tax_levels = taxonomy info for 'entry' (id) in taxonomy_csv dict. (list). i.e. [species, subfamily, ..., order]  # IS THIS NOT EXACTLY THE SAME AS 'TAXONOMY' VARIABLE?
                tax_id = ""                                                # Set tax_id = empty  #IS THIS LINE RELEVANT? SEE COMMENT ON LINE 418

            c_lineage = []  # ADDED
            n = 0

            while tax_id == "" and n < len(tax_levels):

                tax_name = tax_levels[n]
                tax_id = return_ncbi_taxid(entry, tax_name, email_address)
                n += 1
                c_lineage.append(tax_name)

            #lineage_ncbi[entry] = ""
            taxids[entry] = tax_id

            if tax_id == "":                                                   # if tax_id is still empty...
                print(f" - For entry '{entry}' no tax id was found.")
                lineage_custom[entry] = "; ".join(tax_levels)      # the value for the 'entry' key in lineage_custom dict = "tax_levels[0]; tax_levels[1]; ..."

            else:                                                             # if tax_id is not empty...
                #ncbi_l = return_ncbi_lineage(tax_id, email_address)           # ncbi_1 = id obtained from NCBI search for tax_id.
                #lineage_ncbi[entry] = ncbi_l                                  # value for the 'entry' key in lineage_ncbi dict = id obtained from NCBI search.
                lineage_custom[entry] = "; ".join(c_lineage[0:-1])                             # value for the 'entry' key in lineage_custom dict = c_lineage  #? ISN'T THIS EMPTY?

        combined_lineage[entry] = [taxids[entry], lineage_custom[entry]]   # Add to combined_lineage dict key value pair: (entry: [[taxid info],
#      # Need to go through function noting what has been added to each dictionary in each case.
    return combined_lineage


def rejecting_entries(ncbi_lineage, genbank_dict, csv_df, rejection):
    """Print rejected entries to CSV and GenBank files.

    Return df and gb_dict with accepted entries only.
    """
    # ncbi_lineage, genbank_dict, csv_df, rejection = [lineages, new_gb_dict, df_new_ids, args.reject_custom_lineage]

    rejected = []
    new_entries = []
    csv_df.set_index("name", inplace=True)
    # csv_df.rename(columns = {"order" : "taxon"}, inplace=True)          # Rename the "order" column to "taxon". Modify the DataFrame in place

    for record, ncbi_info in ncbi_lineage.items():  # for each key/value pair in ncbi_lineage...

        if ncbi_info[1] != "" and rejection == "True":
            # If the entry has some custom lineage specifications and user wants to reject it...
            rejected.append(record)
            entry = (csv_df.loc[record]).values.tolist()  # Set entry = list of the entries of the "key" row of csv_df
            entry.insert(0, record)  # insert key as first element of entry list
            new_entries.append(entry)  # add entry list to new_entries list

            csv_df.drop([record], inplace=True)  # Drop record (key) column/row (?) from csv_df.

    if len(rejected):
        # Create new DataFrame of rejected entries and drop them from returned DataFrame
        print("\nPrinting rejected entries to CSV file...")

        new_dataframe = pd.DataFrame(new_entries, columns=['name', 'db_id', 'specimen', 'morphospecies', 'species', 'subfamily', 'family', 'order', 'taxid', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'size', 'habitat', 'feeding_behaviour', 'locomotion', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes'])
        del new_dataframe['db_id']
        new_dataframe.to_csv('rejected_metadata.csv', index=False)  # Write new_dataframe to a comma-separated values (csv) file.

        for x in rejected:
            print(f" - Entry '{str(x)}' added to CSV file.")

        # Create GenBank file of rejected entries and drop them from returned genbank_dict
        print("\nPrinting rejected entries to GenBank file...")

        rejected_gb_list = [genbank_dict[x] for x in rejected]
        rejected_gb = open("rejected_entries.gb", 'w')
        SeqIO.write(rejected_gb_list, rejected_gb, "genbank")
        rejected_gb.close()

        for x in rejected:
            print(f" - Entry '{str(x)}' added to GenBank file.")
            del genbank_dict[x]

        print("\nTo keep entries with custom lineage information, run rejected entries again with flag '-r False'.")

    else:
        print("No entries rejected. To reject entries with custom lineage information, re-run script with flag '-r True'.")

    return genbank_dict, csv_df


def insert_taxid(ncbi_lineage, genbank_dict):
    """Insert tax id into gb data (returned from "ncbi_taxid").
    """
    # ncbi_lineage, genbank_dict = [lineages, dict_accepted]

    gb_taxonomy = taxonomy_from_gb(genbank_dict)             # gb_taxonomy = {record_1.name: tax_id, ...}

    for record, tax_id in gb_taxonomy.items():
        genbank_record = genbank_dict[record]                # let genbank_record = the gb info for record.
        ncbi_info = ncbi_lineage[record]                     # let ncbi_info = ncbi lineage info for that record.
        ncbi_id = ncbi_info[0]                               # let ncbi_id = the first element of ncbi_info (the ncbi taxid).
        # -> If there is a taxid in gb - replace it, otherwise create new field for taxid from ncbi or delete if no tax_id given

        if tax_id == "":
            #If no tax_id is given in gb file but is in ncbi_lineage, NCBI taxid will be inserted
            if ncbi_id != "":
                field_given = 0

                for (index, feature) in enumerate(genbank_record.features):                # for all the features of the record...

                    if feature.type == "source":                                           # if the feature type = "source"...
                        #if there is already a "source" field, add db_xref to qualifiers
                        feature.qualifiers['db_xref'] = ['taxon:' + ncbi_id]                      # Set the 'db_xref' qualifier of this feature to 'taxon: first_element_in_ncbi_linage_info'
                        field_given = 1                                                           # Set field_given = 1

                if field_given == 0:
                    #if there is no "source" field, add one and add db_xref to qualifiers
                    len_record = len(genbank_record.seq)                                      # then let len_record = record's sequence length
                    feature_location = FeatureLocation(0, len_record)                         # let feature_location = the region of sequence from 0-END OF SEQ?
                    new_feature = SeqFeature(feature_location, type = 'source')               # let new_feature = make seq
                    new_feature.qualifiers['db_xref'] = ['taxon:' + ncbi_id]
                    genbank_record.features.append(new_feature)

        else:                                                                  # if tax_id is in gb_taxonomy...
            #if there is a taxid in the genbank file
            for (index, feature) in enumerate(genbank_record.features):            # for the features of genbank_record...

                if feature.type == "source":                                           # if the feature type is 'source'...
                    if ncbi_id != "":                                                      # then if ncbi_id is not empty...
                        #replace the gb taxid with the ncbi taxid if there is one
                        feature.qualifiers["db_xref"] = ['taxon:' + ncbi_id]               # repla the 'db_xref' qualifier of this feature to 'taxon: first_element_in_ncbi_linage_info'
                    else:                                                              # if the feature type is not 'source'...
                        #otherwise delete the db_xref qualifier
                        del feature.qualifiers["db_xref"]                                  # delete the 'db_xref qualifier

    return genbank_dict


def loadnamevariants(source=None):
    """Generates dict of name variants for each gene.
    """
    variants = {}
    types = {}
    products = {}

    # Identify source
    if (source is None):
        url = 'https://raw.githubusercontent.com/tjcreedy/biotools/master/gene_name_variants.txt'
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

    # Close handle
    source.close()

    return variants, types, products


def alter_features(genbank_dict):
    """Edit the features in the genbank entries.
    """
    unidentifiable_features = set()
    variants, types, products = loadnamevariants()

    for gb_record, record in genbank_dict.items():                       # for key/value pairs in genbank_dict...

        for (index, feature) in enumerate(record.features):                   # for each of the records features...

            if feature.type.upper() == "CDS":                                        # if the feature's type is "CDS" or "cds"...   ## WAS: if feature.type == "CDS" or feature.type == "cds":
                keys = feature.qualifiers.keys()                                         # then set keys = keys of CDS features.qualifiers dict
                del_features = []                                                        # del_features = []

                for key in keys:                                                         # for each key of CDS qualifiers dict...
                    if key not in ["gene", "location", "codon_start", "transl_table", "label", "product"]:   # if key is not "gene"/"location"/"codon_start"/"trnsl_table"/"label"/"product"...
                        #WHY JUST THESE^? WHY NOT 'translation', 'protein_id' etc.?
                        del_features.append(key)                                                                # add key to del_features list (so this is now a list of all keys in features.qualifiers dict that aren't mentioned above)

                for f in del_features:                                               # for each of these non-mentioned keys...
                    del feature.qualifiers[f]                                            # delete them from features.qualifiers dict

                nametags = ['gene', 'product', 'label', 'standard_name']
                #HOW COULD THERE POSSIBLY BE 'standard_name' WHEN IT WAS DROPPED IN THE LAST STEP?

                if any(t in feature.qualifiers.keys() for t in nametags):               # if there are any CDS qualifier keys left that are also in nametags...
                    name = 0
                    for t in nametags:                                                   # then for those t's that are in nametags and also qualifier keys
                        if t in feature.qualifiers.keys():
                            name = feature.qualifiers[t][0].upper()                         # set name = the value for that qualifier key in uppercase. e.g. for key="gene", name = "NAD6"
                            break

                    if name in variants.keys():                                   # if it (name) is a key in the different_names dict

                        new_name = variants[name]                                     # set new_name = the value for that key    e.g. if name=nad1 or ND1, new_name = ND1

                        feature.qualifiers["gene"] = [new_name]                              # then set it as the value for "gene" in  feature.qualifiers dict    e.g. "gene": "ND1", ..
                        feature.qualifiers["label"] = [f"{new_name} {types[new_name]}"]       # then set the second element of its value in default_qualifier_names as the value for "label" in  feature.qualifiers dict   e.g. "label": "ND1 CDS"
                        feature.qualifiers["product"] = [products[new_name]]      # then set the third element of its value in default_qualifier_names as the value for "product" in  feature.qualifiers dict   e.g. "product": "NADH dehydrogenase subunit 2"

                    else:                                                                # if name is not a key in the different_names dict
                        sys.exit(f"ERROR: Unknown gene name for '{str(gb_record)}' in CDS features: '{str(name)}'")

                else:
                    unidentifiable_features.add((feature.type, feature.location.start, feature.location.end))

        if len(unidentifiable_features):
            sys.stderr.write("\nWARNING: The following sequence entries had unidentifiable annotations:\n")
            for unidfeats in unidentifiable_features:
                sys.stderr.write(gb_record + ": " + ', '.join([f + " " + str(s) + "-" + str(e) for f, s, e in unidfeats]) + "\n")

    return genbank_dict


def add_lineage_df(csv_dataframe, combined_lineage):
    """Add columns with tax_id, custom_ and ncbi_lineage to metadata dataframe.
    """
    # csv_dataframe, combined_lineage = [df_accepted, lineages]

    df_add = pd.DataFrame.from_dict(combined_lineage, orient='index')           # write combined_lineage dict into dataframe called 'df_add' with keys as the index
    df_add.columns = ["taxon_id", "custom_lineage"]                  # change column headers to "taxid", "ncbi_lineage", and "custom_lineage"
    csv_dataframe.drop(['species', 'subfamily', 'family', 'order', 'taxid'], axis=1, inplace=True)                                                    # delete "taxid" column/row from input csv_dataframe
    df = pd.merge(df_add, csv_dataframe, left_index=True, right_index=True)   # merge 'df_add' with 'csv_dataframe', using the index from df_add (keys of combined_lineage dict) and csv_dataframe as the join key(s), calling the resulting dataframe 'df'.
    df.reset_index(level=0, inplace=True)                                     # Remove the level 0 from the index, modifying the dataframe in place.
    df.rename(columns={"index": "name"}, inplace=True)                       # Rename "index" column to "name", modifying the dataframe in place.

    return df

def reformat_df_cols(df):

    #df.rename(columns={'name': 'old_name', 'db_id': 'name'}, inplace=True)

    df = df[['name', 'db_id', 'morphospecies', 'taxon_id', 'custom_lineage', 'specimen', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'authors', 'size', 'habitat', 'feeding_behaviour', 'locomotion', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes']]

    return df


def load_gb_dict_into_db(genbank_data):
    """Load genbank_data as a dictionary into the mysql database.

    db_driver = "MySQLdb"
    db_passwd = "mmgdatabase"
    db_host = "localhost"
    db_user = "root"
    db_name = "mmg_test"
    mysql_engine = "mysql+mysqldb://root:mmgdatabase@localhost/mmg_test"
    namespace = "mmg"
    """
    #genbank_data = records

    print("\nLoading genbank entries into the database...")

    server = BioSeqDatabase.open_database(driver = db_driver, user = db_user, passwd = db_passwd, host = db_host, db = db_name)   # driver = "MySQLdb", user = "root", passwd = "mmgdatabase", host = "localhost", db = "mmg_test"
    db =  server[namespace]
    count = db.load(genbank_data.values())
    server.commit()                             #Commit to memory/save

    print(" - %i sequences loaded." % count)

    return()


def load_df_into_db(csv_dataframe):
    """Loading pandas dataframe with metadata into the database.

    mysql_engine = "mysql+mysqldb://root:mmgdatabase@localhost/mmg_test"
    """
    #csv_dataframe = gb_df_new_ids

    print("\nLoading metadata into the database...")

    engine = create_engine(mysql_engine, echo = False)       # Create an engine object based on URL: "mysql+mysqldb://root:mmgdatabase@localhost/mmg_test", NOT logging the statements as well as a repr() of their parameter lists to the default log handler.
    csv_dataframe.to_sql(name='metadata', if_exists='append', index=False, con=engine)    # write csv_dataframe to an sql database called 'metadata', inserting values to the existing table if it already exists, NOT writing DataFrame index as a column. CON??

    print(" - %i entries loaded." % len(csv_dataframe.index))

    return()



#--------------------


def sql_cols(table, cols, spec):

    #cols = ['taxon_name.name', 'length', 'node_rank']
    #cols = '*'
    #table = 'biosequence'
    #spec = None
    #table, cols, spec = [None, ['count'], ['country=United Kingdom', 'length<25000']]
    #table, cols, spec = [None, ['name', 'db_id'], ['species=cucujiformia', 'length<25000', 'country=United Kingdom']]
    #table, cols, spec = [None, ['*'], ['species=Stenus boops', 'length<25000', 'country!=United Kingdom']]
    #table, cols, spec = [None, ['count'], ['species=Stenus boops']]
    # 'length<25000', 'country!=United Kingdom']]

    # table, cols, spec = [None, None, ['length>25000', 'country=United Kingdom', 'locomotion=arboreal', 'size=12mm']]

    #Reformat inputs
    if spec is None:
        spec = []
    if cols is None:
        cols = []

    #spec = [f"{re.split('=|>|<', s)[0]}{re.findall('=|>|<', s)[0]}{re.split('=|>|<', s)[1]}" if re.split('=|>|<', s)[1].isnumeric() else f"{re.split('=|>|<', s)[0]}='{re.split('=|>|<', s)[1]}'" for s in spec]
    #cols = list(cols)
    all_cols = list(set(cols + [re.split('=|!=|>|<', s)[0] for s in spec]))

    # Unique cols of each table (shared cols assigned to a prioritised table)
    metadata_cols = ['name', 'db_id', 'morphospecies', 'taxon_id', 'custom_lineage', 'specimen', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'metadata.authors', 'size', 'habitat', 'feeding_behavior', 'locomotion', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes']
    bioentry_cols = ['bioentry_id', 'bioentry.biodatabase_id', 'bioentry.taxon_id', 'bioentry.name', 'bioentry.accession', 'identifier', 'division', 'description', 'version']
    bioentry_dbxref_cols = ['bioentry_dbxref.bioentry_id', 'bioentry_dbxref.dbxref_id', 'bioentry_dbxref.rank']
    bioentry_qualifier_value_cols = ['bioentry_qualifier_value.bioentry_id', 'bioentry_qualifier_value.term_id', 'value', 'bioentry_qualifier_value.rank']
    bioentry_reference_cols = ['bioentry_reference.bioentry_id', 'bioentry_reference.reference_id', 'bioentry_reference.start_pos', 'bioentry_reference.end_pos', 'bioentry_reference.rank']
    biosequence_cols = ['biosequence.version', 'length', 'alphabet', 'seq']
    comment_cols = ['comment_id', 'comment.bioentry_id', 'comment_text', 'comment.rank']
    seqfeature_cols = ['seqfeature_id', 'seqfeature.bioentry_id', 'type_term_id', 'source_term_id', 'display_name', 'seqfeature.rank']
    taxon_cols = ['taxon.taxon_id', 'ncbi_taxon_id', 'parent_taxon_id', 'node_rank', 'genetic_code', 'mito_genetic_code', 'left_value', 'right_value']

    #Special queries
    taxonomy = ['subspecies', 'species', 'genus', 'tribe', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']

    #Construct columns dictionary (adding prefixes for table joins)
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

    #Construct tables list
    tables = list(filter(None, list(set([x.split('.')[0] for x in cols_dict.values()] + [table]))))

    #Construct columns string
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
    """Consructs MySQL command
    """
    n = 0
    table_string = start
    while n < len(table_list):
        new_join = f" JOIN {table_list[n]} ON {main_table}.{shared_col}={table_list[n]}.{shared_col}"
        table_string += new_join
        n += 1

    return table_string


def sql_table(tables):

    if len(tables) == 1:

        table_string = tables[0]

    elif len(tables) > 1:

        #Lists of tables sharing linking columns (For duplicates it must be decided which table the column shouls be assigned to).
        BIOENTRY_ID = ['bioentry', 'bioentry_dbxref', 'bioentry_qualifier_value', 'bioentry_reference', 'biosequence', 'comment', 'seqfeature']
        TAXON_ID = ['taxon', 'taxon_name']

        #Split provided tables into groups according to shared columns
        bios = list(set(tables) & set(BIOENTRY_ID))
        taxons = list(set(tables) & set(TAXON_ID))

        ##JOIN TABLES
        joins = ["metadata"]

        #Bios
        if len(bios) >= 1:
            if 'bioentry' in bios:
                bios.remove('bioentry')
            bios_join = table_join(" JOIN bioentry ON metadata.db_id=bioentry.name", bios, 'bioentry', 'bioentry_id')
            joins.append(bios_join)

        #Taxons
        if len(taxons) >= 1:
            if 'taxon' in taxons:
                taxons.remove('taxon')
            taxons_join = table_join(" JOIN taxon ON metadata.taxon_id=taxon.ncbi_taxon_id", taxons, 'taxon', 'taxon_id')
            joins.append(taxons_join)

        table_string = ''.join(joins)

        #COMPLETE_TABLE_STRING = 'metadata JOIN bioentry ON metadata.db_id=bioentry.name JOIN bioentry_dbxref ON bioentry.bioentry_id=bioentry_dbxref.bioentry_id JOIN bioentry_qualifier_value ON bioentry_qualifier_value.bioentry_id=bioentry.bioentry_id JOIN bioentry_reference ON bioentry_reference.bioentry_id=bioentry.bioentry_id JOIN biosequence ON biosequence.bioentry_id=bioentry.bioentry_id JOIN comment ON comment.bioentry_id=bioentry.bioentry_id JOIN seqfeature ON seqfeature.bioentry_id=bioentry.bioentry_id JOIN seqfeature_dbxref ON seqfeature_dbxref.seqfeature_id=seqfeature.seqfeature_id JOIN seqfeature_qualifier_value ON seqfeature_qualifier_value.seqfeature_id=seqfeature.seqfeature_id JOIN taxon ON metadata.taxon_id=taxon.ncbi_taxon_id JOIN taxon_name ON taxon.taxon_id=taxon_name.taxon_id' \

    else:
        sys.exit("ERROR: Cannot construct table. Invalid information provided.")

    return table_string


def sql_spec(tables, cols_dict, spec, spec_type):
    # spec = ['country!=United Kingdom', 'description=Lucanus sp. BMNH 1425267 mitochondrion, complete genome']
    # spec_type='output'
    spec = [s if re.split('=|!=|>|<', s)[1].isnumeric() else f"{re.split('=|!=|>|<', s)[0]}{re.findall('=|!=|>|<', s)[0]}'{re.split('=|!=|>|<', s)[1]}'" for s in spec]

    if len(spec) == 0:
        spec = ''
    else:
        specs = []
        for x in spec:
            s = re.split('=|!=|>|<', x)
            if s[0] in ['subspecies', 'species', 'genus', 'tribe', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']:
                specs.append(f"metadata.taxon_id IN (SELECT DISTINCT include.ncbi_taxon_id FROM taxon INNER JOIN taxon AS include ON (include.left_value BETWEEN taxon.left_value AND taxon.right_value) WHERE taxon.taxon_id IN (SELECT taxon_id FROM taxon_name WHERE name COLLATE LATIN1_GENERAL_CI LIKE '%{s[1][1:-1]}%'))")
            else:
                specs.append(re.findall('=|!=|>|<', x)[0].join([cols_dict[s[0]], s[1]]))

        if spec_type == "output":
            spec = f" WHERE ({') AND ('.join(specs)})"

        if spec_type == "update":
            spec = ', '.join(specs)

    return spec


def construct_sql_output_query(table, cols, spec):
    #table, cols, spec

    tables, cols_string, cols_dict, spec = sql_cols(table, cols, spec)

    table_string = sql_table(tables)

    spec_string = sql_spec(tables, cols_dict, spec, "output")

    mysql_command = f"SELECT {cols_string} FROM {table_string}{spec_string};"

    return mysql_command


def construct_sql_update_query(table, update, spec):

    #update = ['locomotion=arboreal', 'size=12mm']
    #spec = ['length>25000', 'country=United Kingdom']

    tables, _, cols_dict, _ = sql_cols(None, None, spec + update)

    table_string = sql_table(tables)

    spec_string = sql_spec(tables, cols_dict, spec, "output")

    update_string = sql_spec(tables, cols_dict, update, "update")

    mysql_command = f"UPDATE {table_string} SET {update_string}{spec_string};"

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


def fetch_recs(names_dict, db_un, db_pw):
    """Fetches a list of SeqRecords from an input dict of record names/db ids
    """
    # names_dict = {'MH404113': 'GB001', 'KT876913': 'GB007', 'KF364622': 'GB008', 'KT876903': 'GB014'}
    recs = {}

    server = BioSeqDatabase.open_database(driver=db_driver, user=db_un, passwd=db_pw, host=db_host, db=db_name)  # driver = "MySQLdb", user = "root", passwd = "mmgdatabase", host = "localhost", db = "mmg_test"
    db = server[namespace]

    for name, db_id in names_dict.items():
        recs[name] = db.lookup(name=db_id)

    #server.close()

    return recs


def execute_query(mysql_query, db_un, db_pw):
    """Connect to db and execute mysql query
    """
    con = mdb.connect(host=db_host, user=db_un, passwd=db_pw, db=db_name)

    with con:
        cur = con.cursor()
        cur.execute(mysql_query)

    return

"""
def extract_gene(recs, gene_name):
    Extract specific genes from records

    #gene_name = 'nd1'

    #LOOK THROUGH MANUAL ON OTHER WAYS TO DO THIS (ISN'T THERE A BIOSQL AUTO ALTERNATIVE?)
    # SOME RECORDS HAVE NO GENE QUALIFIERS

    subrecs = {}

    #Slice down record
    for idd, record in recs.items():
        for feature in record.features:
            if feature.type.upper() == "CDS" and 'gene' in feature.qualifiers and feature.qualifiers['gene'][0].upper() == gene_name.upper():
                subrec = feature.location.extract(record)
                subrec.description = re.sub(',|complete genome', '', record.description)
                subrecs[idd] = subrec

    #Check for missing genes
    if len(recs) > len(subrecs):
        print(f"WARNING: There are records satisfying your specification that do not contain the '{gene_name.upper()}' gene.")

    return subrecs


def extract_CDS(recs):

    genes = {}

    for x in ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']:
        genes[x] = extract_gene(recs, x)

    return genes
"""

def extract_genes(recs, genes):
    """Extracts genes from SeqRecord objects and writes to dict: {gene_name :  list_of_sliced_seqrecords, ...}
    """
    if '*' in genes:
        genes = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

    print(f"Extracting {len(genes)} genes: {', '.join(genes)}")

    subrecs = {}

    for gene in genes:
        extracted_genes = []
        for idd, record in recs.items():
            for feature in record.features:
                if feature.type.upper() == "CDS" and 'gene' in feature.qualifiers and feature.qualifiers['gene'][0].upper() == gene.upper():  # .upper() not relevant, as alter_features means all in db would be upper, and choice in argparse means all in genes will be upper
                    subrec = feature.location.extract(record)
                    subrec.description = re.sub(', [a-z]+ genome$', '', record.description)
                    subrec.id = record.id
                    subrec.name = record.name
                    extracted_genes.append(subrec)
        subrecs[gene] = extracted_genes

    return subrecs


def csv_from_sql(mysql_command, csv_name, db_un, db_pw):
    #cols, tablename, csv_name, specs = [None, "metadata", "metadata_output", "subregion='Sabah', collectionmethod='BEATING'"]

    #Connect to database and execute command
    con = mdb.connect(host="localhost", user=db_un, passwd=db_pw, db=db_name)
    cur = con.cursor()
    cur.execute(mysql_command)

    #Write data to CSV file
    with open(f"{csv_name}.csv", "w", newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        headers = [i[0] for i in cur.description if i[0] != 'metadata_id']
        csv_writer.writerow(headers)  #Write headers
        csv_writer.writerows(x[1:] for x in cur)

    cur.close()

    return


def seqfile_from_sql(recs_dict, file_name, frmat):
    #Writes list of SeqRecords to a file of chosen format

    #Specific genes
    if any(isinstance(x, list) for x in recs_dict.values()):
        for gene in recs_dict.keys():
            SeqIO.write(recs_dict[gene], f"{file_name}_{gene}.{frmat}", frmat)

    #Full genome
    else:
        SeqIO.write(recs_dict.values(), f"{file_name}.{frmat}", frmat)

    return


def return_count(mysql_command, db_un, db_pw):

    #mysql_command = "SELECT COUNT(*) FROM biosequence WHERE (length<25000);"

    con = mdb.connect(host="localhost", user=db_un, passwd=db_pw, db=db_name)
    cur = con.cursor()
    cur.execute(mysql_command)

    for row in cur:
        print(row[0])

    cur.close()

    return


def mysql_replace_into(table, conn, keys, data_iter):
    """Mimics the behavior of sql REPLACE INTO by passing this callable to to_sql
    """
    from sqlalchemy.dialects.mysql import insert
    from sqlalchemy.ext.compiler import compiles
    from sqlalchemy.sql.expression import Insert
    @compiles(Insert)
    def replace_string(insert, compiler, **kw):
        s = compiler.visit_insert(insert, **kw)
        s = s.replace("INSERT INTO", "REPLACE INTO")
        return s
    data = [dict(zip(keys, row)) for row in data_iter]
    conn.execute(table.table.insert(replace_string=""), data)


def mysql_replace_into_2(table, conn, keys, data_iter):
    """Mimics the behavior of INSERT ... ON DUPLICATE KEY UPDATE ...
    """
    from sqlalchemy.dialects.mysql import insert
    data = [dict(zip(keys, row)) for row in data_iter]
    stmt = insert(table.table).values(data)
    update_stmt = stmt.on_duplicate_key_update(**dict(zip(stmt.inserted.keys(), stmt.inserted.values())))
    conn.execute(update_stmt)


def to_sql_update(df, engine, schema, table):
    """grab a table from a MySQL database (whole or filtered), update/add some rows, and want to perform a REPLACE INTO statement with df.to_sql().
    """
    #df.reset_index(inplace=True)
    sql = f"SELECT column_name FROM information_schema.columns WHERE table_schema = '{schema}' AND table_name = '{table}' AND COLUMN_KEY = 'PRI';"
    id_cols = [x[0] for x in engine.execute(sql).fetchall()]
    id_vals = [df[col_name].tolist() for col_name in id_cols]
    sql = f" DELETE FROM {schema}.{table} WHERE 0"
    for row in zip(*id_vals):
        sql_row = ' AND '.join([" {}='{}' ".format(n, v) for n, v in zip(id_cols, row)])
        sql += f' OR ({sql_row}) '
    engine.execute(sql)
    #df.to_sql(df, engine, schema=schema, if_exists='append', index=False)
    df.to_sql(name=table, if_exists='append', schema=schema, index=False, con=engine)
    return


def overwrite_data(metadata, gb_dict):
    """Overwrite records in the database
    """
    #2 cases:
    #1) Something pulled from db (under db_id), edited, and reingested.
    #2) A new version of a record (under name) needs to be pushed in

    #Grabbing db_id with sql alchemy now futile

    con = mdb.connect(host="localhost", user=db_user, passwd=db_passwd, db=db_name)
    cur = con.cursor()
    server = BioSeqDatabase.open_database(driver=db_driver, user=db_user, passwd=db_passwd, host=db_host, db=db_name)  # driver = "MySQLdb", user = "root", passwd = "mmgdatabase", host = "localhost", db = "mmg_test"
    db = server[namespace]
    engine = create_engine(mysql_engine, echo=False)

    if gb_dict:

        for nom in gb_dict.keys():
            cur.execute(f"SELECT db_id FROM metadata WHERE name='{nom}';")
            for row in cur:
                db_id = row[0]

            for record in db.values():
                if record.name == db_id:
                    idd = record.id
                    accession = record.annotations["accessions"]
                    sequence_version = record.annotations["sequence_version"][0]
                    record = gb_dict[nom]
                    record.name = db_id
                    record.id = idd
                    record.annotations["accessions"] = accession
                    record.annotations["sequence_version"][0] = int(sequence_version) + 1

    if metadata is not None:
        #https://stackoverflow.com/questions/34661318/replace-rows-in-mysql-database-table-with-pandas-dataframe
        #https://tedboy.github.io/pandas/_modules/pandas/io/sql.html

        for nom in metadata['name']:
            cur.execute(f"SELECT db_id FROM metadata WHERE name='{nom}';")
            for row in cur:
                db_id = row[0]

            #Grabs row from local df corresponding to name
            row = metadata.loc[metadata['name'] == nom]
            to_sql_update(row, engine, db_name, 'metadata')
            cur.execute(f"UPDATE metadata SET db_id='{db_id}' WHERE name='{nom}';")

            #BUG: REPLACES DB_ID.

    cur.close()
    server.commit()
    server.close()

    return




# OR row.to_sql(name='metadata', if_exists='append', method=mysql_replace_into_2, index=False, con=engine)





"""

sql_row = pd.io.sql.SQLTable(row, engine, index=False)

row.to_sql(name='metadata', if_exists='append', method=mysql_replace_into_2(sql_row, cur, list(row.columns), row_vals), con=engine)
row = csv_df.loc[csv_df['name'] == ['BIOD00197', 'BIOD01431', 'CCCP00081']]
row = csv_df.loc[['BIOD00197', 'BIOD01431', 'CCCP00081']]

row_vals = row.values.tolist()
row_vals = [['QINL077', None, 'MG73', None, 'Galerucinae', 'Chrysomelidae', 'Coleoptera', None, None, 'Adult', None, None, 'Qinling', 'China', np.nan, np.nan, ';)', None, None, None, 'Nie et al', 'lib7', None, 'Qinling Chrysomelidae', None, None, None]]

engine = create_engine(mysql_engine, echo=False)

row.to_sql(name='metadata', if_exists='append', method=mysql_replace_into_2(row, cur, list(row.columns), row_vals), con=engine)

row.to_sql(name='metadata', if_exists='append', method=mysql_replace_into_2, index=False, con=engine)

sql_row = pd.io.sql.SQLTable(row, engine, index=False)


#For mapping Pandas tables to SQL tables.
'BIOD02014'

sql_row = pd.io.sql.SQLTable(row)


help(pd.io.sql.SQLTable)
"""


def mysql_replace_into(table, conn, keys, data_iter):
    #Mimics the behavior of sql REPLACE INTO by passing this callable to to_sql
    
    from sqlalchemy.dialects.mysql import insert
    from sqlalchemy.ext.compiler import compiles
    from sqlalchemy.sql.expression import Insert
    @compiles(Insert)
    def replace_string(insert, compiler, **kw):
        s = compiler.visit_insert(insert, **kw)
        s = s.replace("INSERT INTO", "REPLACE INTO")
        return s
    data = [dict(zip(keys, row)) for row in data_iter]
    conn.execute(table.insert(replace_string=""), data)



def mysql_replace_into_2(table, conn, keys, data_iter):
    #Mimics the behavior of INSERT ... ON DUPLICATE KEY UPDATE ...
    
    #keys = cols
    #data_iter = row_vals

    from sqlalchemy.dialects.mysql import insert
    data = [dict(zip(keys, row)) for row in data_iter]
    stmt = insert(table.table).values(data)
    update_stmt = stmt.on_duplicate_key_update(**dict(zip(stmt.inserted.keys(), stmt.inserted.values())))
    conn.execute(update_stmt)
