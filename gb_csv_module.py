#!/urs/bin/env python

"""Functions to modify and upload a genbank-file and metadata in a csv-file into a MySQL database.
"""


###IMPORT PACKAGES

import sys, time, urllib.request

##Import modules for handeling the genbank files:
from Bio import SeqIO, Entrez, SeqFeature                             # Entrez is a molecular biology database system that provides integrated access to NCBI’s databases such as PubMed, GenBank, GEO, and many others.
from Bio.SeqFeature import SeqFeature, FeatureLocation

##Import modules for handeling the csv files:
import pandas as pd
import numpy as np

##Import modules for dealing with the mysql database:
import MySQLdb as mdb
from BioSQL import BioSeqDatabase
from sqlalchemy import create_engine


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


def correct_header(csv_dataframe):
    """Check metadata file has correct columns #by creating a list of its column headers and checking it against a list of the expected headers.
    """
    # csv_dataframe = new_csv_df

    csv_header = csv_dataframe.columns.values.tolist()                                     # LUKE: creates list of csv column headers.      ## Pandas index.tolist() Converts dataframes column values to list, like below.  ###6 LUKE: So why use the .values function when .columns returns an index already?
    expected_header = ['name', 'specimen', 'morphospecies', 'species', 'subfamily', 'family', 'order', 'taxid', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes']

    if expected_header != csv_header:
        print("Incorrect header in CSV file.\n")
        sys.exit("Current header is: " + str(csv_header) + "\n\nIt must be as follows: ['name', 'specimen', 'morphospecies', 'species', 'subfamily', 'family', 'order', 'taxid', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes']")

    return()


def matching_inputids(csv_dataframe, gb_dictionary):
    """Check if the GenBank and CSV metadata files have matching entries, and that there are no duplicates in the CSV file. (Duplicates in the GenBank file would
    have generated an error when first loaded into a dictionary in gb_into_dictionary.)
    """
    # gb_dictionary, csv_dataframe = [gb_dict, csv_df]

    new_gb_dict, new_csv_df = [gb_dictionary, csv_dataframe]

    ids_csv = csv_dataframe['name'].values.tolist()
    unique = list(set(ids_csv + list(gb_dictionary.keys())))                        #all IDs in both files (the union).
    shared = list(set(ids_csv) & set(gb_dictionary.keys()))                         #all IDs shared by by both files (the intersection).
    discrepant_ids = [a for a in unique if a not in shared]                         #all IDs that aren't shared by both files.
    csv_duplicates = set([idd for idd in ids_csv if ids_csv.count(idd) > 1])  #all duplicate IDs in the CSV file.

    if len(csv_duplicates):
        # Check IDs are unique in the CSV metadata file.
        sys.exit("ERROR: There are multiple rows sharing the same name in your CSV file: " + ', '.join(csv_duplicates) + ". IDs must be unique.")

    if len(shared) != len(unique):
        # If the IDS in the CSV and GenBank files are not identical...
        x = input("Your CSV and GenBank files contain different entries.\nWould you like to ignore these and proceed with shared entries only ('P') or cancel the operation ('C')?\n?>").capitalize()

        while not (x == 'C' or x == 'P'):
            x = input("Type 'P' to ignore discrepant entries and proceed, or type 'C' to cancel the operation.\n?>").capitalize()

        if x == 'C':
            sys.exit("Operation cancelled.")

        else:
            # return new gb_dict with discrepant entries deleted
            for gb_record in gb_dictionary:
                if gb_record in discrepant_ids:
                    print(f" - Skipping entry '{gb_record}' as it appears in the GenBank file but not the CSV file.")

            new_gb_dict = {key: gb_dictionary[key] for key in gb_dictionary if key not in discrepant_ids}

            # return new csv_dataframe with discrepant entries deleted
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

    if len(str(startvalue)) > int(padding):
        sys.exit("ERROR: The starting number " + str(startvalue) + " exceeds the digits for 0-padding (" + str(padding) + " digits).")

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


def check_new_ids(dict_new_ids):
    """Check if the new database ids already exist in the database.
    """

    for old_id, new_id in dict_new_ids.items():           #for each key (old id) in the dict_new_ids dictionary...
        # new_id = the value for that key. (Are the keys the old ids and the values the new ones?)
        mysql_command = "SELECT * FROM bioentry WHERE name = '" + str(new_id) + "'"     # set mysql_command = the command: Select everything from the bioentry database where the name = new_id
        con = mdb.connect(host = db_host, user = db_user, passwd = db_passwd, db = db_name)    #Connects to MySQL database

        with con:                          # With the connection as defined above...
            cur = con.cursor()             # Create a cursor (allows Python code to execute SQL commands in a database session) and bind it to the connection for the entire lifetime (this is what con.cursor() does).
            cur.execute(mysql_command)     # Execute the mysql_command command. (.execute() executes a SQL statement.)
            result = cur.fetchone()        # Fetch the next row of the query (i.e. selected) result set, returning a single tuple, or None when no more data is available

        if result is not None:             # If data is still available... (i.e. if there are any entries in the bioentry database with name = new_id...)
            sys.exit(f"The new database id '{str(new_id)}' is already present in the database:\n{str(result)}")

    return()


def change_ids_genbank(genbank_dict, dict_new_ids, key):
    """Change the ids in the genbank file to the new database ids (LLLLNNNNN).
    """

    for gb_record, record in genbank_dict.items():                                                                    # L: for each key in genbank_dict, let 'record' equal the value for that key.

        #Change the name (record.name) (to key) of each entry to new database id (to value)
        record.name = dict_new_ids[record.name]                                                                     # record_name = record.name = dict_new_ids[record_name] = new_name

        #Change the accession (record.id) to the new database id including version number (.0)
        new_accession = record.name + ".0"
        record.id = new_accession
        record.annotations["accessions"] = [record.name]                                                 #.annotations is a Python dictionary.

        #Change the version
        record.annotations["sequence_version"] = 0       # Genbank files include a version annotation. This sets it to 0, because we're going to use that ourselves to record our own versions of an entry.

        if key == "DESCRIPTION":
            #Erase the description as it contained the input id
            record.description = ""

    return()


def change_names_csv(csv_dataframe, dict_new_ids):
    """Create dictionary with new database ids and old input names to df.
    """

    dict_df = pd.DataFrame(list(dict_new_ids.items()), columns = ["name", "db_id"])  # list(dict_new_ids.items()) creates a list of tuples of the dict pairs - (key, value)
    #creates a 2-column dataframe with the "name" column featuring all the the old ids and the "db_id" column featuring the new ids.

    #Merge csv_dataframe and dict_df -> append new database ids as a column
    new_csv_df = pd.merge(dict_df, csv_dataframe, on = 'name')                       # this merges both dataframes on the 'name' (old id) column, which in effect simply adds the "db_id" (new id) column in dict_df onto the csv_dataframe

    return new_csv_df  # returns csv_dataframe with the new ids added on in a single new column called "db_id"

def change_names_gb_csv(csv_dataframe, dict_new_ids):
    """Create dictionary with new database ids and old input names to df.
    """

    dict_df = pd.DataFrame(list(dict_new_ids.items()), columns = ["accession", "db_id"])  # list(dict_new_ids.items()) creates a list of tuples of the dict pairs - (key, value)
    #creates a 2-column dataframe with the "name" column featuring all the the old ids and the "db_id" column featuring the new ids.

    #Merge csv_dataframe and dict_df -> append new database ids as a column
    new_csv_df = pd.merge(dict_df, csv_dataframe, on = 'accession')                       # this merges both dataframes on the 'name' (old id) column, which in effect simply adds the "db_id" (new id) column in dict_df onto the csv_dataframe

    return new_csv_df

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


def taxonomy_from_gb(genbank_dict):
    """Function returning dictionary of taxids from genbank file if it has it (or "" if not).
    """

    gb_taxa = {}

    for gb_record, record in genbank_dict.items():

        tax_id = ""

        for (index, feature) in enumerate(record.features):                                                                                # Need to understand more about record attributes and their functions for this. (Features are SeqFeature objects)

            if feature.type == "source":        # .type - the specified type of the feature (ie. CDS, exon, repeat...)

                if "db_xref" in feature.qualifiers:       # THIS SHOULD INSTEADC BE 'DBXREFS'???      qualifiers - A dictionary of qualifiers on the feature (holds metadata about a feature if it is present in the entry). These are analogous to the qualifiers from a GenBank feature table. The keys of the dictionary are qualifier names, the values are the qualifier values.
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

        new_dataframe = pd.DataFrame(new_entries, columns=['name', 'db_id', 'specimen', 'morphospecies', 'species', 'subfamily', 'family', 'order', 'taxid', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes'])
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


def loadnamevariants():
    # generates dict of name variants for each gene.
    output = {}
    for line in urllib.request.urlopen("https://raw.githubusercontent.com/tjcreedy/biotools/master/gene_name_variants.txt"):
        line = line.decode('utf-8').strip()
        name = line.split(";")[0]
        annotype = line.split(":")[0].split(";")[1]
        variants = line.split(":")[1].split(",")
        for v in variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' ' + annotype.upper()]:
                    output[v + s] = name
    return output


def alter_features(genbank_dict):
    # Edit the features in the genbank entries.

    unidentifiable_features = set()
    different_names = loadnamevariants()
    default_qualifier_names = {"ND1": ["ND1", "ND1 CDS", "NADH dehydrogenase subunit 1"],
                               "ND2": ["ND2", "ND2 CDS", "NADH dehydrogenase subunit 2"],
                               "ND3": ["ND3", "ND3 CDS", "NADH dehydrogenase subunit 3"],
                               "ND4": ["ND4"," ND4 CDS", "NADH dehydrogenase subunit 4"],
                               "ND4L": ["ND4L", "ND4L CDS", "NADH dehydrogenase subunit 4L"],
                               "ND5": ["ND5", "ND5 CDS", "NADH dehydrogenase subunit 5"],
                               "ND6": ["ND6", "ND6 CDS", "NADH dehydrogenase subunit 6"],
                               "CYTB": ["CYTB", "CYTB CDS", "cytochrome b"],
                               "COX1": ["COX1", "COX1 CDS", "cytochrome c oxdiase subunit 1"],
                               "COX2": ["COX2", "COX2 CDS", "cytochrome c oxdiase subunit 2"],
                               "COX3": ["COX3", "COX3 CDS", "cytochrome c oxdiase subunit 3"],
                               "ATP6": ["ATP6", "ATP6 CDS", "ATP synthase F0 subunit 6"],
                               "ATP8": ["ATP8", "ATP8 CDS", "ATP synthase F0 subunit 8"]}

    for gb_record, record in genbank_dict.items():                       # for key/value pairs in genbank_dict...

        for (index, feature) in enumerate(record.features):                   # for each of the records features...

            if feature.type.upper() == "CDS":                                        # if the feature's type is "CDS" or "cds"...   ## WAS: if feature.type == "CDS" or feature.type == "cds":
                keys = feature.qualifiers.keys()                                         # then set keys = keys of CDS features.qualifiers dict
                del_features = []                                                        # del_features = []

                for key in keys:                                                         # for each key of CDS qualifiers dict...
                    if key not in ["gene", "location", "codon_start", "trnsl_table", "label", "product"]:   # if key is not "gene"/"location"/"codon_start"/"trnsl_table"/"label"/"product"...
                        del_features.append(key)                                                                # add key to del_features list (so this is now a list of all keys in features.qualifiers dict that aren't mentioned above)

                for f in del_features:                                               # for each of these non-mentioned keys...
                    del feature.qualifiers[f]                                            # delete them from features.qualifiers dict

                nametags = ['gene', 'product', 'label', 'standard_name']

                if any(t in feature.qualifiers.keys() for t in nametags):               # if there are any CDS qualifier keys left that are also in nametags...
                    name = 0
                    for t in nametags:                                                   # then for those t's that are in nametags and also qualifier keys
                        if t in feature.qualifiers.keys():
                            name = feature.qualifiers[t][0].upper()                         # set name = the value for that qualifier key in uppercase. e.g. for key="gene", name = "NAD6"
                            break

                    if name in different_names.keys():                                   # if it (name) is a key in the different_names dict

                        new_name = different_names[name]                                     # set new_name = the value for that key    e.g. if name=nad1 or ND1, new_name = ND1

                        feature.qualifiers["gene"] = new_name                                # then set it as the value for "gene" in  feature.qualifiers dict    e.g. "gene": "ND1", ..
                        feature.qualifiers["label"] = default_qualifier_names[new_name][1]        # then set the second element of its value in default_qualifier_names as the value for "label" in  feature.qualifiers dict   e.g. "label": "ND1 CDS"
                        feature.qualifiers["product"] = default_qualifier_names[new_name][2]      # then set the third element of its value in default_qualifier_names as the value for "product" in  feature.qualifiers dict   e.g. "product": "NADH dehydrogenase subunit 2"

                    else:                                                                # if name is not a key in the different_names dict
                        sys.exit("Unknown gene name for " + str(gb_record) + " in CDS features: " + str(name))

                else:
                    unidentifiable_features.add((feature.type, feature.location.start, feature.location.end))

        if len(unidentifiable_features) > 0:
            sys.stderr.write("\nWARNING\nThe following sequence entries had unidentifiable annotations:\n")
            for unidfeats in unidentifiable_features:
                sys.stderr.write(gb_record + ": " + ', '.join([f + " " + str(s) + "-" + str(e) for f, s, e in unidfeats]) + "\n")

    return genbank_dict


def add_lineage_df(csv_dataframe, combined_lineage):
    """Add columns with tax_id, custom_ and ncbi_lineage to metadata dataframe.
    """
    # csv_dataframe, combined_lineage = [df_accepted, lineages]

    df_add = pd.DataFrame.from_dict(combined_lineage, orient = 'index')           # write combined_lineage dict into dataframe called 'df_add' with keys as the index
    df_add.columns = ["taxon_id", "custom_lineage"]                  # change column headers to "taxid", "ncbi_lineage", and "custom_lineage"
    csv_dataframe.drop(['species', 'subfamily', 'family', 'order', 'taxid'], axis=1, inplace=True)                                                    # delete "taxid" column/row from input csv_dataframe
    df = pd.merge(df_add, csv_dataframe, left_index = True, right_index = True)   # merge 'df_add' with 'csv_dataframe', using the index from df_add (keys of combined_lineage dict) and csv_dataframe as the join key(s), calling the resulting dataframe 'df'.
    df.reset_index(level = 0, inplace = True)                                     # Remove the level 0 from the index, modifying the dataframe in place.
    df.rename(columns = {"index" : "name"}, inplace = True)                       # Rename "index" column to "name", modifying the dataframe in place.

    return df

def reorder_df_cols(df):

    df = df[['name', 'db_id', 'morphospecies', 'taxon_id', 'custom_lineage', 'specimen', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes']]

    return()


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
    db =  server[namespace]                    # db = server["mmg"]       .new_database(namespace, description='Testing MMG')
    count = db.load(genbank_data.values())      # load dict values (ONLY?) into database
    server.commit()                             # Commit to memory/save

    print(" - %i sequences loaded." % count)        # is count an integer?

    return()


def load_df_into_db(csv_dataframe):
    """Loading pandas dataframe with metadata into the database.

    mysql_engine = "mysql+mysqldb://root:mmgdatabase@localhost/mmg_test"
    """
    #csv_dataframe = gb_df_new_ids

    print("\nLoading metadata into the database...")

    engine = create_engine(mysql_engine, echo = False)       # Create an engine object based on URL: "mysql+mysqldb://root:mmgdatabase@localhost/mmg_test", NOT logging the statements as well as a repr() of their parameter lists to the default log handler.
    csv_dataframe.to_sql(name = 'metadata', if_exists = 'append', index = False, con = engine)    # write csv_dataframe to an sql database called 'metadata', inserting values to the existing table if it already exists, NOT writing DataFrame index as a column. CON??

    print(" - %i entries loaded." % len(csv_dataframe.index))

    return()








#---------------------


def text_to_list(accessions):
    """Converts text-file of IDs (one per line) into a list, exiting if duplicates are present.
    """
    #accessions = args.input_accessions

    #Create a comma-delimited list of accession numbers from text file (stripping any blank spaces/empty lines).
    with open(accessions, "r") as acc:
        accs = acc.read()
        striplist = lambda lis:[x.strip() for x in lis]
        acc_list = list(filter(None, striplist(accs.split('\n'))))

    #Check for duplicates
    duplicates = set([idd for idd in acc_list if acc_list.count(idd) > 1])

    if len(duplicates):
        x = input("ERROR: There are duplicate IDs in your text-file: '" + "', '".join(duplicates) + "'. IDs must be unique.\nWould you like to skip these and proceed ('P') or cancel the operation ('C')?\n?>").capitalize()
        while not (x == 'P' or x == 'C'):
            x = input("Type 'P' to proceed with unique IDs only or 'C' to cancel the operation.\n?>").capitalize()
        if x == 'C':
            sys.exit("Operation cancelled.")
        else:
            print("Skipping entries: '" + "', '".join(duplicates) + "'.")
            acc_list = [a for a in acc_list if a not in duplicates]

    return acc_list


def chunker(seq, size):

    if type(seq) is set:
        seq = list(seq)

    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def return_gb_data(acc_list, email):
    """Takes list of of IDs/accessions and returns dict of corresponding GenBank entries.
    """
    #email = args.users_email

    records = {}
    Entrez.email = email

    for accset in chunker(acc_list, 200):
        with Entrez.efetch(db='nuccore', id=accset, rettype='gb') as handle:
            for record in SeqIO.parse(handle, "gb"):
                records[record.name] = record
        time.sleep(0.1)

    return records


def extract_metadata(records):
    """Extracts metadata from gb_dict and writes to DataFrame.
    """

    #Extract metadata to separate dict
    gb_metadata = {}
    for id, record in records.items():
        for feature in record.features:
            if feature.type == "source":
                taxid = str(feature.qualifiers['db_xref'])[8:-2]
                gb_metadata[id] = taxid

    #Write to DataFrame
    gb_met_df = pd.DataFrame.from_dict(gb_metadata, orient='index')
    gb_met_df.reset_index(inplace=True)
    gb_met_df.columns = ["accession", "taxon_id"]

    for label in ['name', 'morphospecies', 'custom_lineage', 'specimen', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'authors', 'library', 'datasubmitter', 'projectname', 'uin', 'notes']:
        gb_met_df[label] = ""

    for header in ['latitude', 'longitude']:
        gb_met_df[header] = np.nan

    return gb_met_df

