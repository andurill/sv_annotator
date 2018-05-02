#!/usr/bin/env python

import os.path
import pandas as pd
import sys


def make_cytoband_table(file):
    '''
    creates a panda dataframe of chromosome cytoband mapping
    '''
    # type: (File) -> pandaDataFrame
    try:
        cb_df = pd.read_csv(file, sep="\t", \
        usecols=["Chr", "Arm", "Band", "Bp_start", "Bp_stop"])
    except IOError as e:
        print e
        raise
    return cb_df


def make_canonicalTx_table(file):
    '''
    creates a dictionary of canonical transcript mapped to Gene name
    '''
    # type: (File) -> Dict[Str, Str]
    Tx_dict = {}
    try:
        with open(file, r) as f:
            for line in f:
                if "#" in line:
                    continue
                else:
                    line = line.strip().split("\t")
                    gene = line[0]
                    tx_id = line[1]
                    Tx_dict[tx_id] = gene
    except IOError as e:
        print e
        raise
    return Tx_dict


def make_refFlatTable(file):
    '''
    creates a panda dataframe of all canonical transcripts
    '''
    # type: (File) -> pandaDataFrame
    try:
        refFlat_df = pd.read_csv(file, sep="\t", header=0, \
        column=["Chr", "Start", "End", "Strand", "Gene", "Trascript",\
        "Exon", "cDNA", "AA"])
    except IOError as e:
        print e
        raise
    return refFlat_df


def make_kinaseGeneTable(file):
    '''
    creates a panda dataframe of all kinase genes
    '''
    # type: (File) -> pandaDataFrame
    try:
        kinase_df = pd.read_csv(file, sep="\t", header=0, \
                                usecols=["Chr", "Start", "End", "Strand",\
                                "Refseq_transcriptID", "HUGO", "Kinase_domain_Type",\
                                "Exon_start", "Exon_end", "AA_start", "AA_end",\
                                "cDNA_start", "cDNA_end"])
    except IOError as e:
        print e
        raise
    return kinase_df