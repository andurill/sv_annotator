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
    creates a pandas data frame of refseq canonical transcripts
    mapped to ensembl transcripts and hugo name
    '''
    # type: (File) -> pandaDataFrame
    try:
        Tx_df = pd.read_csv(file, sep="\t", header=0, \
        column=["Ensembl", "Refseq", "Gene"])
    except IOError as e:
        print e
        raise
    return Tx_df


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