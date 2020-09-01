#!/usr/bin/env python

"""Functions to analyze genbank-files and csv-files."""

## Imports ##
import sys
import pandas as pd
import numpy as np


## Functions ##

##For each entry get sequence length, topology and feature count from gb-file

def gb_length(genbank_dict):
    """Get the sequence lengths for each entry.
    Returns a dict containing key value pairs - ID: seqlength
    """
    lengths = {}

    for gb_record in genbank_dict:
        record = genbank_dict[gb_record]
        name = record.name
        length = len(record.seq)
        lengths[name] = length

    return lengths

def gb_topology(genbank_dict):
    """Get the topology (linear vs circular) for each entry.
    Returns a dict containing key value pairs - ID: topology
    """
    topologies = {}

    for gb_name, gb_record in genbank_dict.items():
        # gb_name, gb_record = list(genbank_dict.items())[0]
        # print("NAME: " + str(gb_name) + "\nREC.NAME: " + str(gb_record.name))
        if 'data_file_division' not in gb_record.annotations.keys():
            sys.exit(
                str(gb_name) + " does not have a topology annotation, wtf?")
        else:
            topology = gb_record.annotations["data_file_division"]
            topologies[gb_record.name] = topology

    return topologies

def gb_features(genbank_dict):
    """Get the features count (gene, cds, trna, rrna) for each entry.
    Returns a dict containing key value pairs.

     Output
        - {name: [gene_count, cds_count, trna_count, rrna_count], ...}
    """
    features = {}

    for gb_record, record in genbank_dict.items():

        gene_count = 0
        cds_count = 0
        trna_count = 0
        rrna_count = 0

        for (index, feature) in enumerate(record.features):

            if feature.type == "gene":
                gene_count += 1
            if feature.type == "CDS":
                cds_count += 1
            if feature.type == "tRNA":
                trna_count += trna_count
            if feature.type == "rRNA":
                rrna_count += rrna_count

        features[record.name] = [gene_count, cds_count, trna_count, rrna_count]

    return features

##Write sequence length, topology and feature count to a dataframe

def dict_to_df(lengths, topologies, features):
    """Combine sequence lengths, topologies and feature counts of entries into
    a pandas dataframe.
    """
    df_l = pd.DataFrame.from_dict(lengths, orient='index')
    df_l.columns = ["length"]

    df_t = pd.DataFrame.from_dict(topologies, orient='index')
    df_t.columns = ["topology"]

    df_f = pd.DataFrame.from_dict(features, orient='index')
    df_f.columns = ["gene", "cds", "trna", "rrna"]

    df = pd.merge(pd.merge(df_l, df_t, left_index=True, right_index=True), df_f,
                  left_index=True, right_index=True)
    df.reset_index(level=0, inplace=True)
    df.rename(columns={"index": "name"}, inplace=True)

    df.to_csv('features.csv', index=False, header=True)

    return df

##Completeness of sequences: Calculate percentages for the genbank featuers,
# seq length and topology

def calculations(df_features):
    """Count and calculate percentages for different features.
    """
    f = open("summary.txt", "w")

    # Sequence length
    l_s = df_features.length
    l_a = l_s.values
    count = (l_a >= 14000).sum()
    smaller = (l_a < 14000).sum()
    larger = (l_a > 17000).sum()
    total = count + smaller
    count_p = round((count / total) * 100, 2)
    smaller_p = round((smaller / total) * 100, 2)
    larger_p = round((larger / total) * 100, 2)
    f.write("Summary\n\nTotal number of sequences: " + str(total) + "\n")
    f.write("Sequence length 14,000 - 17,000 bp:\t" + str(
        count - larger) + "\t(" + str(count_p) + "%)\n")
    f.write("Sequence length smaller than 14,000 bp:\t" + str(
        smaller) + "\t(" + str(smaller_p) + "%)\n")
    f.write(
        "Sequence length larger than 17,000 bp:\t" + str(larger) + "\t(" + str(
            larger_p) + "%)\n\n")

    # Topology (linear vxs circular)
    t_a = df_features.topology
    t_s = t_a.values
    linear = (t_s == "linear").sum()
    linear_p = round((linear / total) * 100, 2)
    circular = (t_s == "circular").sum()
    circular_p = round((circular / total) * 100, 2)
    f.write(
        "Topology:\nEntries with linear genomes:\t" + str(linear) + "\t(" + str(
            linear_p) + "%)\n")
    f.write("Entries with circular genomes:\t" + str(circular) + "\t(" + str(
        circular_p) + "%)\n\n")

    # Gene number
    g_s = df_features.gene
    g_a = g_s.values
    genes = (g_a == 13).sum()
    less = (g_a < 13).sum()
    more = (g_a > 13).sum()
    genes_p = round((genes / total) * 100, 2)
    less_p = round((less / total) * 100, 2)
    more_p = round((more / total) * 100, 2)
    f.write("Gene number:\nNumber of genes is 13:\t" + str(genes) + "\t(" + str(
        genes_p) + "%)\n")
    f.write("More than 13 genes:\t" + str(more) + "\t(" + str(more_p) + "%)\n")
    f.write(
        "Fewer than 13 genes:\t" + str(less) + "\t(" + str(less_p) + "%)\n\n")

    # CDSs number
    c_s = df_features.cds
    c_a = c_s.values
    cds = (c_a == 13).sum()
    lessc = (c_a < 13).sum()
    morec = (c_a > 13).sum()
    cds_p = round((cds / total) * 100, 2)
    lessc_p = round((lessc / total) * 100, 2)
    morec_p = round((morec / total) * 100, 2)
    f.write("CDS number:\nNumber of CDSs is 13:\t" + str(cds) + "\t(" + str(
        cds_p) + "%)\n")
    f.write("More than 13 CDSs:\t" + str(morec) + "\t(" + str(morec_p) + "%)\n")
    f.write(
        "Fewer than 13 CDSs:\t" + str(lessc) + "\t(" + str(lessc_p) + "%)\n\n")

    # tRNA number
    tr_s = df_features.trna
    tr_a = tr_s.values
    trnas = (tr_a == 22).sum()
    lesstr = (tr_a < 22).sum()
    moretr = (tr_a > 22).sum()
    trnas_p = round((trnas / total) * 100, 2)
    lesstr_p = round((lesstr / total) * 100, 2)
    moretr_p = round((moretr / total) * 100, 2)
    f.write("tRNA number:\nNumber of tRNAs is 22:\t" + str(trnas) + "\t(" + str(
        trnas_p) + "%)\n")
    f.write(
        "More than 22 tRNAs:\t" + str(moretr) + "\t(" + str(moretr_p) + "%)\n")
    f.write("Fewer than 22 tRNAs:\t" + str(lesstr) + "\t(" + str(
        lesstr_p) + "%)\n\n")

    # rRNA number
    sr_s = df_features.rrna
    sr_a = sr_s.values
    srnas = (sr_a == 2).sum()
    lesssr = (sr_a < 2).sum()
    moresr = (sr_a > 2).sum()
    srnas_p = round((srnas / total) * 100, 2)
    lesssr_p = round((lesssr / total) * 100, 2)
    moresr_p = round((moresr / total) * 100, 2)
    f.write("rRNA number:\nNumber of rRNAs is 2:\t" + str(srnas) + "\t(" + str(
        srnas_p) + "%)\n")
    f.write(
        "More than 2 rRNAs:\t" + str(moresr) + "\t(" + str(moresr_p) + "%)\n")
    f.write("Fewer than 2 rRNAs:\t" + str(lesssr) + "\t(" + str(
        lesssr_p) + "%)\n\n")

    f.close()

    return

## For seq length and each feature calculate mean, std, median, min and max
# values

def calculations_ftr(df_features):
    """Calculations for different features.
    """
    column = df_features.length
    lengths = column.values
    column = df_features.cds
    cdss = column.values
    column = df_features.gene
    genes = column.values
    column = df_features.trna
    trnas = column.values
    column = df_features.rrna
    rrnas = column.values

    f = open("summary.txt", "a")
    f.write("\nCalculations\n\n")

    feature_names = ["Sequence length", "Genes", "CDSs", "tRNAs", "rRNAs"]

    for ftr, name in zip([lengths, genes, cdss, trnas, rrnas], feature_names):
        m = round(np.nanmean(ftr), 2)
        std = round(np.nanstd(ftr), 2)
        md = round(np.nanmedian(ftr), 2)
        mi = np.nanmin(ftr)
        ma = np.nanmax(ftr)

        f.write(str(name) + ":\n")
        f.write(
            "mean = " + str(m) + ", std = " + str(std) + ", median = " + str(
                md) + ", range = " + str(mi) + " - " + str(ma) + "\n")

    f.close()

    return

##For each gene (CDS) in the genbank entries analyze sequence length, extract
# lengths and output as dataframe

def contig_length(genbank_dict):
    """Create dictionary with a dictionary for each entry containing the
    sequence length of each CDS & look at variation.
    
    Output as pandas dataframe.
    """
    # genbank_dict = gb_dict

    gene_lengths = {}

    for gb_record in genbank_dict:

        record = genbank_dict[gb_record]
        name = record.name
        feature_length = {}

        for (index, feature) in enumerate(record.features):

            if feature.type.upper() == "CDS":

                if "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                    cds_len = len(feature.location.extract(record).seq)
                    feature_length[gene_name] = int(cds_len)

        gene_lengths[name] = feature_length

    df = pd.DataFrame.from_dict(gene_lengths, orient='index')
    df.reset_index(level=0, inplace=True)
    df.rename(columns={"index": "entry"}, inplace=True)
    df.to_csv('contig_lengths.csv', index=False, header=True)

    return df

def length_variation(gene_lengths):
    """Take a dataframe with lengths for each cds.
    """
    gene_names = gene_lengths.columns.values.tolist()
    del gene_names[0]

    f = open("summary.txt", "a")
    f.write("\n\nGene length variation\n\n")

    for gene in gene_names:
        std_m = []
        column = gene_lengths.loc[:, gene]
        lengths = np.array(column)

        m = round(np.nanmean(lengths), 2)
        std = round(np.nanstd(lengths), 2)
        var = round(np.nanvar(lengths), 2)

        md = round(np.nanmedian(lengths), 2)
        mi = np.nanmin(lengths)
        ma = np.nanmax(lengths)

        f.write(str(gene) + ":\n")
        f.write("mean = " + str(m) + ", std = " + str(std) + ", var = " + str(
            var) + ", median = " + str(md) + ", range = " + str(
            mi) + " - " + str(ma) + "\n")

    f.close()

    return

##Analyze gene order

def gene_order(genbank_dict):
    """Extract the start and end position of each sequence for each gene.
    """
    gene_order = {}

    for gb_record in genbank_dict:
        record = genbank_dict[gb_record]
        name = record.name

        genes = {}

        for (index, feature) in enumerate(record.features):
            position = []

            if feature.type == "CDS" or feature.type == "cds":
                gene = feature.qualifiers["gene"][0]
                start = feature.location._start.position
                end = feature.location._end.position
                position.append(start)
                position.append(end)
                genes[gene] = position

        gene_order[name] = genes

    df = pd.DataFrame.from_dict(gene_order, orient='index')
    df.reset_index(level=0, inplace=True)
    df.rename(columns={"index": "entry"}, inplace=True)

    return df

def isexactsubset(query, subject):
    """Check if the query is part of the subject.
    """
    if set(query).issubset(set(subject)):
        issubset = [subject[i:i + len(query)] == query or subject[i:i - len(
            query):-1] == query for i, x in enumerate(subject) if x == query[0]]
        return any(issubset)

    else:
        return False

def extract_order(genbank_dict):
    """Extract the genes in the order they are listed in the genbank features.
    """
    gene_order = {}

    for gb_record in genbank_dict:

        record = genbank_dict[gb_record]
        name = record.name
        genes = []

        for (index, feature) in enumerate(record.features):

            if feature.type == "CDS" or feature.type == "cds":
                gene = feature.qualifiers["gene"]
                genes.append(gene)

        gene_order[name] = tuple(genes)

    return gene_order

def analyze_order(order):
    """Analyze if the gene order is correct and contiguous.
    """
    ancestral_order = ("COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND5",
                       "ND4", "ND4L", "ND6", "CYTB", "ND1", "ND2", "COX1",
                       "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND5", "ND4",
                       "ND4L", "ND6", "CYTB", "ND1", "ND2")
    different_order = {}

    for entry in order:
        genes = order[entry]
        if genes:
            if not isexactsubset(genes, ancestral_order):
                different_order[entry] = genes

    f = open("summary.txt", "a")
    f.write("\n\nGene order\n\n")
    f.write("Different gene order is given for: " + ', '.join(
        list(different_order.keys())) + '.')
    f.close()

    return different_order

##Analyze family and country distribution

def family_distribution(metadata_df):
    """Extract families from metadata table.
    """
    # metadata_df=csv_df
    families = metadata_df.family
    strip = families.str.strip()
    family_count = strip.value_counts()

    family_count.to_csv('family_count.csv', index_label='Family',
                        header=['Count'])

    return

def country_distribution(metadata_df):
    """Extract countries from metadata table.
    """
    # metadata_df = csv_df
    countries = metadata_df.country
    country_count = countries.value_counts(dropna=False)
    country_count.to_csv('country_count.csv', index_label='Country',
                         header=['Count'])

    return ()

##Save gene sequences in one fasta file

def gene_sequence(genbank_dict):
    """Save genes sequences in fasta-file (one file for complete sequences
    only, one for both complete or partial sequences).
    """
    f = open("gene_sequences.fasta", 'w')
    fc = open("gene_sequences.complete_only.fasta", 'w')

    for gb_record, record in genbank_dict.items():

        name = record.name

        for (index, feature) in enumerate(record.features):

            if feature.type == "CDS" or feature.type == "cds":

                if "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"]
                    cds_seq = feature.location.extract(record).seq
                    loc = feature.location
                    if not ">" in str(loc) and not "<" in str(loc):
                        fc.write(">" + str(name) + "_" + str(gene_name) + "\n")
                        fc.write(str(cds_seq) + "\n")

                    f.write(">" + str(name) + "_" + str(gene_name) + "\n")
                    f.write(str(cds_seq) + "\n")

    fc.close()
    f.close()

    return

"""
BUGS

BUG1: gb_topology() - 'data_file_division' and 'topology' mixup - Line 39
      Also affects Line 137-138 in summary.txt
 """
