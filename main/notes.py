import os, sys
import pandas as pd
import numpy as np
from notes_constants import *


refFlat_summary = pd.read_csv("./data/refFlat_summary.txt", sep="\t", dtype={
                              'a': str, 'b': str, 'c': str, 'd': str, 
                              'e': int, 'f': int, 'g': int})


def get_bkp_info(bkp):
    #bkp.firstexon = refFlat_summary[refFlat_summary.Gene ==
    #                               "EGFR"].first_exon.tolist().pop(0)
    bkp_dict = refFlat_summary[
        (refFlat_summary['Gene'].values == bkp.gene) & 
        (refFlat_summary['Transcript'].values == bkp.Transcript)]\
        .to_dict('records').pop()
    bkp.firstexon = bkp_dict['first_exon']
    bkp.lastexon = bkp_dict['last_exon']
    bkp.exon, bkp.intron = [None]*2
    if bkp.desc.starswith("Exon "):
        bkp.exon = int(bkp.desc.split(" ")[1])
    elif bkp.desc.starswith("Intron"):
        exon = bkp.desc.split(" ")[6]
        if "before" in bkp.desc:
            bkp.intron = int(exon) - 1
            bkp.exon = int(exon)
        else:
            bkp.intron = int(exon)
            bkp.exon = int(exon)
    if bkp.intron:
        bkp.site = "intron " + str(bkp.intron)
    else:
        bkp.site = "exon " + str(bkp.exon)
    return


def get_exons_involved(sv):
    if sv.isFusion:
        get_bkp_info(sv.fusionPartner1)
        get_bkp_info(sv.fusionPartner2)
        note1 = "%s exons %s - %s" % (sv.fusionPartner1.gene,
                                      sv.fusionPartner1.firstexon,
                                      sv.fusionPartner1.exon)
        note2 = "%s exons %s - %s" % (sv.fusionPartner2.gene,
                                      sv.fusionPartner2.exon,
                                      sv.fusionPartner2.lastexon)
    elif sv.bkp1.isPanel and sv.bkp2.isPanel and sv.bkp1.isCoding and sv.bkp2.isCoding:
        get_bkp_info(sv.annotationPartner1)
        get_bkp_info(sv.annotationPartner2)
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s %s and %s %s" % (sv.annotationPartner1.gene,
                                         sv.annotationPartner1.site,
                                         sv.annotationPartner2.gene,
                                         sv.annotationPartner2.site)
        else:
            note1 = "%s exons %s - %s" % (sv.annotationPartner1.gene,
                                          sv.annotationPartner1.firstexon,
                                          sv.annotationPartner1.exon)
        note2 = "%s exons %s - %s" % (sv.fusionPartner2.gene,
                                      sv.fusionPartner2.exon,
                                      sv.fusionPartner2.lastexon)
        gene1, tx1, cdna1 = sv.annotationPartner1.gene, \
            sv.annotationPartner1.transcript, sv.annotationPartner1.cdna
        gene2, tx2, cdna2 = sv.annotationPartner2.gene, \
            sv.annotationPartner2.transcript, sv.annotationPartner2.cdna
        if sv.isIntragenic:
            return "%s (%s) %s: %s_%s:%s" %\
                (gene1, tx1, fusion_type, cdna1, cdna2, svtype)
        else:
            return "%s (%s) - %s (%s) %s: %s:%s_%s:%s%s" %\
                (gene1, tx1, gene2, tx2, fusion_type,
                 cdna1, gene1, cdna2, gene2, svtype)
    elif sv.bkp1.isPanel and sv.bkp1.isCoding:
        gene1, tx1, cdna1 = sv.bkp1.gene, sv.bkp1.transcript, sv.bkp1.cdna
        if not cdna1.startswith("chr"):
            cdna2 = "chr" + sv.bkp2.chrom + ":g." + str(sv.bkp2.pos)
            return "%s (%s) %s: %s:%s_%s%s" %\
                (gene1, tx1, fusion_type, cdna1, gene1, cdna2, svtype)
        else:
            raise BreakPointIntergenic(sv.bkp1)
    else:
        gene2, tx2, cdna2 = sv.bkp2.gene, sv.bkp2.transcript, sv.bkp2.cdna
        if not cdna2.startswith("chr"):
            cdna1 = "chr" + sv.bkp1.chrom + ":g." + str(sv.bkp1.pos)
            return "%s (%s) %s: %s:%s_%s%s" %\
                (gene2, tx2, fusion_type, cdna2, gene2, cdna1, svtype)
        else:
            raise BreakPointIntergenic(sv.bkp2)
    return

def get_note(sv):
    note = None
    return note


def get_prefix(sv):
    prefix = "Note: The "
    fusion_type = "rearrangement"
    if sv.isFusion:
        if sv.isKnownFusion:
            fusion_type = "fusion"
        prefix += "%s - %s " % (sv.fusionPartner)
    return


def get_kinase(sv):
    for bkp in (sv.bkp1, sv.bkp2):
        if bkp.isKinase:

    

def functional_significance(sv):
    if sv.isKnownFusion or\
    (sv.isFusion and ((sv.bkp1.isHotspot and sv.bkp1.isEntireKinase) or\
    (sv.bkp2.isHotspot and sv.bkp2.isEntireKinase))) or\
    (sv.svtype == "DELETION" and\
    ((sv.bkp1.isPanel and sv.bkp1.isTumourSuppressor) or\
    (sv.bkp2.isPanel and sv.bkp2.isTumourSuppressor))):
        return ""
    else:
        return " The functional significance is undetermined."