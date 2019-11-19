#!/urs/bin/env python

"""Functions to modify and upload a genbank-file and metadata in a csv-file into a MySQL database.
"""


###IMPORT PACKAGES

import sys, time

##Import modules for handeling the genbank files:
from Bio import SeqIO, Entrez, SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation

##Import modules for handeling the csv files:
import pandas as pd

##Import modules for dealing with the mysql database:
import MySQLdb as mdb
from BioSQL import BioSeqDatabase
from sqlalchemy import create_engine


###VARIABLES

db_driver = "MySQLdb"
db_user = "root"
db_passwd = "mmgdatabase"
db_host = "localhost"
db_name = "mmg_test"
mysql_engine = "mysql+mysqldb://root:mmgdatabase@localhost/mmg_test"

namespace = "mmg"


###DEFINE FUNCTIONS


def gb_into_dictionary(gb_filename, key):
    '''Take a file in genbank-format and load it into a dictionary, define function for keys (default: by name in LOCUS).
    '''
    
    if key == "LOCUS":

        def get_seqname(record):
            #Given a SeqRecord, return the sequence name as a string.
            seqname = str(record.name)
            return(seqname)

        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"), key_function = get_seqname)

    elif key == "ACCESSION":
        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"))

    elif key == "DEFINITION":
        
        def get_definition(record):
            #Given a SeqRecord, return the definition as a string.
            definition = str(record.description)
            return(definition)
        
        gb_dictionary = SeqIO.to_dict(SeqIO.parse(gb_filename, "genbank"), key_function = get_definition)

    return(gb_dictionary)


def correct_header(csv_dataframe):
    '''Check metadata file has correct columns.
    '''
    
    csv_header = csv_dataframe.columns.values.tolist()
    expected_header = ['name' ,'specimen' ,'morphospecies', 'species', 'subfamily', 'family', 'order', 'taxid', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes']
    
    if expected_header != csv_header:
        print("Your header is: " + str(csv_header))
        sys.exit('The header in the csv-file is wrong.\nIt must be as follows: name, specimen, morphospecies, species, subfamily, family, order, taxid, collectionmethod, lifestage, site, locality, subregion, country, latitude, longitude, authors, library, datasubmitter, projectname, accession, uin, notes')
    
    return()


def matching_inputids(csv_dataframe, gb_dictionary):
    '''Check if the genbank file and metadata have matching set of entries by input id.
    '''
    
    ids_csv=csv_dataframe['name'].values.tolist()
    
    if len(ids_csv) != len(gb_dictionary):
        sys.exit("The number of entries in the gengbank- and csv-file does not match.")

    for contigname in ids_csv:
        
        if not contigname in gb_dictionary:
            sys.exit("For entry '" + str(contigname) + "' listed in the csv-file there is no entry in the genbank-file.")

    for gb_record in gb_dictionary:
        
        record_id = gb_dictionary[gb_record].name
        
        if record_id not in ids_csv:
            sys.exit("For entry '" + str(record_id) + "' listed in the gebank-file there is no entry in the csv-file.")

    return()


def unique_csv_ids(csv_dataframe):
    '''Check ids are unique in the metadata file.
    '''
    
    ids_csv = csv_dataframe['name'].values.tolist()
    set_ids_csv = set(ids_csv)
    
    if int(len(ids_csv)) != int(len(set_ids_csv)):
        sys.exit('Duplicated entries in csv-file! IDs must be unique.')
        
    return()


def new_ids(genbank_dict,prefix,integer,padding):
    '''Check if the new ids are looking fine given the input (prefix, integer, padding) by the user.
    '''
    
    dict_new_ids = {}
    
    if len(prefix) > 4:
        sys.exit("The prefix is too long (not more than 4 characters).")

    if prefix.isalpha() == False:
        sys.exit("The prefix should only consist of letters.")

    if str(integer).isdigit() == False:
        sys.exit("The number should only consist of digits.")

    if len(str(len(genbank_dict.keys()))) > int(padding):
        sys.exit("A 0-padding by " + str(padding) + " digits is not sufficient for the ingestion of " + str(len(genbank_dict)) + " new entries.")

    if len(str(integer)) > int(padding):
        sys.exit("The starting number " + str(integer) + " exceeds the digits for 0-padding (" + str(padding) + " digits).")

    for gb_record in genbank_dict:

        record = genbank_dict[gb_record]
        record_name = record.name
        pad = "{0:0" + str(padding) + "d}"
        padded = pad.format(int(integer))
        new_name = prefix + str(padded)
        dict_new_ids[record_name] = new_name
        integer = int(integer) + 1

    return(dict_new_ids)


def check_new_ids(dict_new_ids):
    '''Check if the new database ids already exist in the database.
    '''

    for old_id in dict_new_ids:

        new_id = dict_new_ids[old_id]
        mysql_command = "SELECT * FROM bioentry WHERE name = '" + str(new_id) + "'"
        con = mdb.connect(host = db_host, user = db_user, passwd = db_passwd, db = db_name)

        with con:
            cur = con.cursor()
            cur.execute(mysql_command)
            result = cur.fetchone()

        if result is not None:
            sys.exit("The new database id '" + str(new_id) + "' is already present in the database:\n" + str(result))

    return()


def change_ids_genbank(genbank_dict, dict_new_ids, key):
    '''Change the ids in the genbank file to the new database ids (LLLLNNNNN).
    '''

    for gb_record in genbank_dict:
        record = genbank_dict[gb_record]

        #Change the name (record.name) (to key) of each entry to new database id (to value)
        record_name = record.name
        new_name = dict_new_ids[record_name]
        record.name = new_name

        #Change the accession (record.id) to the new database id including version number (.0)
        new_accession = new_name + ".0"
        record.id = new_accession
        record.annotations["accessions"] = [new_name]

        #Change the version
        record.annotations["sequence_version"] = 0

        if key == "DESCRIPTION":
            #Erase the description as it contained the input id
            record.description = ""

    return()


def change_names_csv(csv_dataframe, dict_new_ids):
    '''Create dictionary with new database ids and old input names to df.
    '''

    dict_df = pd.DataFrame(list(dict_new_ids.items()), columns = ["name", "db_id"])
    
    #Merge csv_dataframe and dict_df -> append new database ids as a column
    new_csv_df = pd.merge(dict_df, csv_dataframe, on = 'name')
    
    return(new_csv_df)


def taxonomy_metadata(csv_dataframe):
    '''Function returning a dictionary with all the taxonomic information for each id from metadata ([] if no info given).
    '''
    
    csv_taxa = {}
    ids_csv = csv_dataframe['name'].values.tolist()
    df = csv_dataframe.set_index("name")

    for id in ids_csv:

        tax_list = []
        species = df.loc[id, 'species']
        subfamily = df.loc[id, 'subfamily']
        family = df.loc[id, 'family']
        order = df.loc[id, 'order']

        if not pd.isnull(species):
            tax_list.append(species)
        if not pd.isnull(subfamily):
            tax_list.append(subfamily)
        if not pd.isnull(family):
            tax_list.append(family)
        if not pd.isnull(order):
            tax_list.append(order)

        csv_taxa[id] = tax_list

    return(csv_taxa)


def taxonomy_from_gb(genbank_dict):
    '''Function returning dictionary of taxids from genbank file if it has it (or "" if not).
    '''
    
    gb_taxa = {}

    for gb_record in genbank_dict:

        record = genbank_dict[gb_record]
        record_name = genbank_dict[gb_record].name
        tax_id = ""

        for (index, feature) in enumerate(record.features):

            if feature.type == "source":

                if "db_xref" in feature.qualifiers:
                    xref = feature.qualifiers["db_xref"]
                    
                    if len(xref) > 1:
                        print("There are several ids in the db_xref 'source' information.")
                    else:
                        tax_id = xref[0].split(":")[1]

        gb_taxa[record_name] = tax_id

    return(gb_taxa)


def taxid_metadata(csv_dataframe):
    '''Return dictionary with taxids provided in metadata table ("" if not given).
    '''
    
    csv_taxids = {}
    ids_csv = csv_dataframe['name'].values.tolist()
    df = csv_dataframe.set_index("name")

    for id in ids_csv:

        taxid = df.loc[id, 'taxid']

        if not pd.isnull(taxid):
            csv_taxids[id] = taxid

        else:
            csv_taxids[id] = ""

    return(csv_taxids)


def return_ncbi_taxid(searchterm, email_address):
    '''For each entry get tax_id from NCBI taxonomy based on taxonomic information.
    '''
    
    Entrez.email = email_address
    handle = Entrez.esearch(db = "taxonomy", retmax = 2, term = searchterm)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]
    
    if len(id_list) == 0:
        print("No hits found for search term '" + str(searchterm) + "'in NCBI taxonomy.")
        tax_id = ""
        
    elif len(id_list) > 1:
        print("Multiple hits found for search term '" + str(searchterm) + "' in NCBI taxonomy.")
        tax_id = ""
   
    elif len(id_list) == 1:
        tax_id = id_list[0]
   
    time.sleep(0.5)

    return(tax_id)


def return_ncbi_lineage(searchterm, email_address):
    '''Search NCBI for lineage informatin given a tax id.
    '''
    
    Entrez.email = email_address

    handle = Entrez.efetch(db = "taxonomy", id = searchterm)
    record = Entrez.read(handle)
    handle.close()

    if len(record) == 0:
        print("No hits found for tax id '" + str(searchterm) + "'in NCBI taxonomy.")

    elif len(record) > 1:
        print("Multiple hits found for tax id '" + str(searchterm) + "' in NCBI taxonomy.")

    elif len(record) == 1:
        taxonomy = record[0]["Lineage"]
        taxon = record[0]["ScientificName"] #Add the taxon itself to the lineage
        lineage = taxonomy + "; " + taxon

    time.sleep(0.5)

    return(lineage)


def ncbi_lineage(csv_dataframe, email_address, searchterm):
    '''Search on NCBI for tax ids if not given and if tax ids given in the first place or found search for lineage.

    Incorporates functions "taxonomy_metadata", "taxid_metadata" & "return_ncbi_taxid".
    '''

    taxonomy_csv = taxonomy_metadata(csv_dataframe) #Get all tax info from csv as dict
    taxids_csv = taxid_metadata(csv_dataframe) #Get all tax ids from csv as dict

    print("Searching NCBI for taxonomy ...")
    combined_lineage = {}
    lineage_ncbi = {}
    lineage_custom = {}
    taxids = {}

    for entry in taxids_csv:
        given_taxid = taxids_csv[entry]

        if given_taxid != "":
            #Search given tax id on ncbi & add ncbi_lineage
            ncbi_l = return_ncbi_lineage(given_taxid, email_address)
            lineage_ncbi[entry] = ncbi_l
            lineage_custom[entry] = ""
            taxids[entry] = given_taxid

        else:
            #No tax id in metadata table
            taxonomy = taxonomy_csv[entry] #this returns a list of tax levels or an empty list []

            if taxonomy == []:
                #No tax info is provided in the csv -> reject or use user input!
                if isinstance(searchterm, str):
                    tax_id = return_ncbi_taxid(searchterm, email_address)
                    tax_levels = [searchterm]

                else:
                    print("For entry '" + str(entry) + "' no information about the taxonomy is given in the csv-file.")
                    term_to_search = input("What term should be searched in NCBI taxonomy?\n")
                    tax_levels = [term_to_search]
                    tax_id = return_ncbi_taxid(term_to_search, email_address)
                    ##If nothing is found the tax id is ""
            else:
                #Tax info is searched on ncbi
                tax_levels = taxonomy_csv[entry]
                tax_id = ""

            custom = ""
            c_lineage = ""
            n = 0

            while tax_id == "" and n < len(tax_levels):

                c_lineage = c_lineage + custom
                tax_name = tax_levels[n]
                tax_id = return_ncbi_taxid(tax_name, email_address)
                n = n + 1
                custom = tax_name + "; "

            lineage_ncbi[entry] = ""
            taxids[entry] = tax_id

            if tax_id == "":
                print("For entry '" + str(entry) + "' no tax id was found.")
                seperator = ";"
                lineage_custom[entry] = seperator.join(tax_levels)

            else:
                ncbi_l = return_ncbi_lineage(tax_id, email_address)
                lineage_ncbi[entry] = ncbi_l
                lineage_custom[entry] = c_lineage

        combined_lineage[entry] = [taxids[entry], lineage_ncbi[entry], lineage_custom[entry]]

    return(combined_lineage)


def rejecting_entries_new_gb(ncbi_lineage, genbank_dict, csv_dataframe, rejection):
    '''Print rejected entries to genbank file.

    Return genbank dict with accepted entries only.
    '''

    rejected_gb = open("rejected_entries.gb", 'w')
    gb_taxonomy = taxonomy_from_gb(genbank_dict)

    for record in gb_taxonomy:
        gb_id = gb_taxonomy[record]
        ncbi_info = ncbi_lineage[record]

        if ncbi_info[2] != "":
            #The entry has some custom lineage specification
            if rejection == "True":
                #The user wants to reject the entry
                SeqIO.write(genbank_dict[record], rejected_gb, "genbank")
                del genbank_dict[record]

    rejected_gb.close()

    return(genbank_dict)


def rejecting_entries_new_csv(ncbi_lineage, csv_df, rejection):
    '''Print rejected entries to csv file.

    Return df with accepted entries only.
    '''

    new_entries = []
    csv_df.set_index("name", inplace=True)
    csv_df.rename(columns = {"order" : "taxon"}, inplace=True)

    for record in ncbi_lineage:
        ncbi_info = ncbi_lineage[record]

        if ncbi_info[2] != "":
            #The entry has some custom lineage specifications
            print("The entry " + str(record) + " has custom lineage information id is provided.")

            if rejection == "True":
                print("To keep entries with custom lineage information, run rejected entries again with flag '-r False'.")
                entry = (csv_df.loc[record]).values.tolist()
                entry.insert(0, record)
                new_entries.append(entry)

                csv_df.drop([record], inplace = True)

    new_dataframe = pd.DataFrame(new_entries, columns = ['name', 'db_id', 'specimen', 'morphospecies', 'species', 'subfamily', 'family', 'order', 'taxid', 'collectionmethod', 'lifestage', 'site', 'locality', 'subregion', 'country', 'latitude', 'longitude', 'authors', 'library', 'datasubmitter', 'projectname', 'accession', 'uin', 'notes'])
    del new_dataframe['db_id']
    new_dataframe.to_csv('rejected_metadata.csv', index = False)

    return(csv_df)


def insert_taxid(ncbi_lineage, genbank_dict):
    '''Insert tax id into gb data (returned from "ncbi_taxid").
    '''

    gb_taxonomy = taxonomy_from_gb(genbank_dict)

    for record in gb_taxonomy:
        genbank_record = genbank_dict[record]
        tax_id = gb_taxonomy[record]
        ncbi_info = ncbi_lineage[record]
        ncbi_id = ncbi_info[0]
        # -> If there is a taxid in gb - replace it, otherwise create new field for taxid from ncbi or delete if no tax_id given

        if tax_id == "":
            #No, gb id given -> NCBI taxid will be inserted if there is one
            if ncbi_id != "":	
                field_given = 0

                for (index, feature) in enumerate(genbank_record.features):

                    if feature.type == "source":
                        feature.qualifiers['db_xref'] = ['taxon:' + ncbi_id]
                        field_given = 1

                if field_given == 0:
                    len_record = len(genbank_record.seq)
                    feature_location = FeatureLocation(0, len_record)
                    new_feature = SeqFeature(feature_location, type = 'source')
                    new_feature.qualifiers['db_xref'] = ['taxon:' + ncbi_id]
                    genbank_record.features.append(new_feature)
                #if there is already a source field or not... add source or only db_xref

        else:
            #There is a taxid in the genbank file
            for (index, feature) in enumerate(genbank_record.features):

                if feature.type == "source":
                    if ncbi_id != "":
                        feature.qualifiers["db_xref"] = ['taxon:' + ncbi_id]
                    else:
                        del feature.qualifiers["db_xref"]

    return(genbank_dict)


def alter_features(genbank_dict):
    '''Edit the features in the genbank entries.
    '''
    
    different_names={"ND1":"ND1", "nad1":"ND1", "ND2":"ND2", "nad2":"ND2", "ND3":"ND3", "nad3":"ND3", "ND4":"ND4", "nad4":"ND4", "ND4L":"ND4L", "nad4l":"ND4L", "ND4l":"ND4L", "ND5":"ND5", "nad5":"ND5", "ND6":"ND6", "nad6":"ND6", "CYTB":"CYTB", "cytb":"CYTB", "COX1":"COX1", "cox1":"COX1", "COI":"COX1", "COX2":"COX2", "cox2":"COX2", "COII":"COX2", "COX3":"COX3", "cox3":"COX3", "COIII":"COX3", "ATP6":"ATP6", "atp6":"ATP6", "ATP8":"ATP8", "atp8":"ATP8"}
    different_features={"ND1":["ND1", "ND1 CDS", "NADH dehydrogenase subunit 1"],"ND2":["ND2", "ND2 CDS", "NADH dehydrogenase subunit 2"],"ND3":["ND3", "ND3 CDS", "NADH dehydrogenase subunit 3"],"ND4":["ND4"," ND4 CDS", "NADH dehydrogenase subunit 4"],"ND4L":["ND4L", "ND4L CDS", "NADH dehydrogenase subunit 4L"],"ND5":["ND5", "ND5 CDS", "NADH dehydrogenase subunit 5"],"ND6":["ND6", "ND6 CDS", "NADH dehydrogenase subunit 6"],"CYTB":["CYTB", "CYTB CDS", "cytochrome b"],"COX1":["COX1", "COX1 CDS", "cytochrome c oxdiase subunit 1"],"COX2":["COX2", "COX2 CDS", "cytochrome c oxdiase subunit 2"],"COX3":["COX3", "COX3 CDS", "cytochrome c oxdiase subunit 3"],"ATP6":["ATP6", "ATP6 CDS", "ATP synthase F0 subunit 6"],"ATP8":["ATP8", "ATP8 CDS", "ATP synthase F0 subunit 8"]}

    for gb_record in genbank_dict:
        record = genbank_dict[gb_record]

        for (index,feature) in enumerate(record.features):

            if feature.type == "CDS" or feature.type == "cds":
                keys = feature.qualifiers.keys()
                del_features = []

                for key in keys:
                    if key != "gene" and key != "location" and key != "codon_start" and key != "trnsl_table" and key != "label" and key != "product":
                        del_features.append(key)

                for f in del_features:
                    del feature.qualifiers[f]

                if "gene" in feature.qualifiers:
                    name = feature.qualifiers["gene"]

                    if type(name) == list:
                        name = feature.qualifiers["gene"][0].split(";")[0]

                    if name in different_names.keys():

                        new_name = different_names[name]

                        feature.qualifiers["gene"] = new_name
                        feature.qualifiers["label"] = different_features[new_name][1]
                        feature.qualifiers["product"] = different_features[new_name][2]

                    else:
                        sys.exit("Unknown gene name for " + str(gb_record) + " in CDS features: " + str(name))

    return(genbank_dict)


def add_lineage_df(csv_dataframe, combined_lineage):
    '''Add columns with tax_id, custom_ and ncbi_lineage to metadata dataframe.
    '''

    df_add = pd.DataFrame.from_dict(combined_lineage, orient = 'index')
    df_add.columns = ["taxid", "ncbi_lineage", "custom_lineage"]
    del csv_dataframe["taxid"]
    df = pd.merge(df_add, csv_dataframe, left_index = True, right_index = True)
    df.reset_index(level = 0, inplace = True)
    df.rename(columns = {"index" : "name"}, inplace = True)

    return(df)


def load_gb_dict_into_db(genbank_data):
    '''Load genbank_data as a dictionary into the mysql database.
    '''
    
    print("Loading genbank entries into the database ...")

    server = BioSeqDatabase.open_database(driver = db_driver, user = db_user, passwd = db_passwd, host = db_host, db = db_name)
    db = server[namespace]
    count = db.load(genbank_data.values())
    server.commit()

    print("Loaded %i sequences" % count)

    return()


def load_df_into_db(csv_dataframe):
    '''Loading pandas dataframe with metadata into the database.
    '''

    print("Loading metadata into the database ...")

    engine = create_engine(mysql_engine, echo = False)
    csv_dataframe.to_sql(name = 'metadata', if_exists = 'append', index = False, con = engine)

    return()

