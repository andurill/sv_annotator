#!usr/bin/env python2

import os
import sys
import pandas as pd
import numpy as np
from notes_constants import *


refFlat_summary = pd.read_csv("./data/refFlat_summary.txt", sep="\t", dtype={
                              'a': str, 'b': str, 'c': str, 'd': str,
                              'e': int, 'f': int, 'g': int})


def get_bkp_info(bkp):
    """
    Get exon and intron features for a breakpoint object.
    bkp -> None
    """
    # bkp.firstexon = refFlat_summary[refFlat_summary.Gene ==
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
    """
    Get exons involved in an sv object based on the variant type
    and the breakpoint sites
    sv -> None
    """
    if sv.isFusion:
        get_bkp_info(sv.fusionPartner1)
        get_bkp_info(sv.fusionPartner2)
        note1 = "%s exons %s - %s" % (sv.fusionPartner1.gene,
                                      sv.fusionPartner1.firstexon,
                                      sv.fusionPartner1.exon)
        note2 = "%s exons %s - %s" % (sv.fusionPartner2.gene,
                                      sv.fusionPartner2.exon,
                                      sv.fusionPartner2.lastexon)
        sv.fusionPartner1.variantSite1 = sv.fusionPartner1.firstexon
        sv.fusionPartner1.variantSite2 = sv.fusionPartner1.exon
        sv.fusionPartner2.variantSite1 = sv.fusionPartner2.firstexon
        sv.fusionPartner2.variantSite2 = sv.fusionPartner2.exon
    elif sv.bkp1.isPanel and sv.bkp2.isPanel and sv.bkp1.isCoding and sv.bkp2.isCoding:
        get_bkp_info(sv.annotationPartner1)
        get_bkp_info(sv.annotationPartner2)
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s %s and %s %s" % (sv.annotationPartner1.gene,
                                         sv.annotationPartner1.site,
                                         sv.annotationPartner2.gene,
                                         sv.annotationPartner2.site)
        elif sv.isIntragenic:
            note1 = "%s exons %s - %s" % (sv.annotationPartner1.gene,
                                          sv.annotationPartner1.exon,
                                          sv.annotationPartner2.exon)
        else:
            note1 = "%s exons %s - %s" % (sv.annotationPartner1.gene,
                                          sv.annotationPartner1.firstexon,
                                          sv.annotationPartner1.exon)
            note2 = "%s exons %s - %s" % (sv.annotationPartner2.gene,
                                          sv.annotationPartner2.exon,
                                          sv.annotationPartner2.lastexon)
    elif sv.bkp1.isPanel and sv.bkp1.isCoding:
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s %s" % (sv.annotationPartner1.gene,
                               sv.annotationPartner1.site)
        elif sv.bkp1.strand == "+":
            note1 = "%s exons %s - %s" % (sv.annotationPartner1.gene,
                                          sv.annotationPartner1.exon,
                                          sv.annotationPartner1.lastexon)
        elif sv.bkp1.strand == "-":
            note1 = "%s exons %s - %s" % (sv.annotationPartner1.gene,
                                          sv.annotationPartner1.firstexon,
                                          sv.annotationPartner1.exon)
    else:
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s %s" % (sv.annotationPartner2.gene,
                               sv.annotationPartner2.site)
        elif sv.bkp2.strand == "+":
            note1 = "%s exons %s - %s" % (sv.annotationPartner2.gene,
                                          sv.annotationPartner2.exon,
                                          sv.annotationPartner2.lastexon)
        elif sv.bkp2.strand == "-":
            note1 = "%s exons %s - %s" % (sv.annotationPartner2.gene,
                                          sv.annotationPartner2.firstexon,
                                          sv.annotationPartner2.exon)
    return


def get_note(sv):
    note = None
    return note


def get_prefix(sv):
    """
    Get the prefix of a note based on variant type
    sv -> str
    """
    prefix = "Note: The "
    fusion_type = "rearrangement"
    if sv.isKnownFusion:
        fusion_type = "fusion"
        prefix += str(sv.annotation.split(":")[0]) + \
            fusion_type + " involves "
    elif sv.isFusion:
        prefix += str(sv.annotation.split(":")[0]) + \
            fusion_type + " is a " + sv.svtype.lower() + \
            "that results in a fusion of "
    elif sv.isIntragenic:
        prefix += str(sv.annotation.split(":")[0]) + \
            fusion_type + " is an intragenic " + \
            sv.svtype.lower() + " of "
    else:
        prefix += str(sv.annotation.split(":")[0]) + \
            fusion_type + " is a " + sv.svtype.lower() + " of "
    return prefix


def get_fusion_misc(sv):
    """
    get frame and kinase domain information for
    fusion variants only
    sv -> str
    """
    for bkp in (sv.bkp1, sv.bkp2):
        if bkp.isKinase:
            return


def functional_significance(sv):
    """
    Determine if functional significace note is neccessary
    based on variant type
    sv -> str
    """
    if sv.isKnownFusion or\
        (sv.isFusion and
         ((sv.fusionPartner1.isHotspot and sv.fusionPartner1.isEntireKinase) or
          (sv.fusionPartner2.isHotspot and sv.fusionPartner2.isEntireKinase))) or\
        (sv.svtype == "DELETION" and
         ((sv.bkp1.isPanel and sv.bkp1.isTumourSuppressor) or
          (sv.bkp2.isPanel and sv.bkp2.isTumourSuppressor))):
        return ""
    else:
        return " The functional significance is undetermined."
