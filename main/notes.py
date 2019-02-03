#!usr/bin/env python2

import os
import sys
import pandas as pd
import numpy as np
from notes_constants import *


refFlat_summary = pd.read_csv(
    "./data/refFlat_summary.txt", sep="\t", 
        dtype={'a': str, 'b': str, 'c': str, 
        'd': str, 'e': int, 'f': int, 'g': int})


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


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def get_kinase_annotation(bkp, kinase_annotation):
    kinase_interval = []
    try:
        interval = (bkp.variantSite1, bkp.variantSite2)
        kinase_interval = kinase_annotation[
            (kinase_annotation['HUGO'].values == bkp.gene)].\
            iloc[0, [1, 4]].values.tolist()
    except IndexError:
            # no values found
    except NameError:
            # not fusion
    except ValueError:
            # what
    if kinase_interval:
        overlap = getOverlap(interal, kinase_interval)
        if overlap == 0:
            return "does not include the kinase domain of %s" % (bkp.gene)
        else:
            if min(interval) <= min(kinase_interval) and \
                max(interval) >= max(kinase_interval):
                return "includes the kinase domain of %s" % (bkp.gene)
            else:
                return "includes a part of the kinase domain of %s" % (bkp.gene)
    return None
        

def get_note(sv):
    note = None
    return note


#" and ".join(filter(lambda x: isinstance(x, str), [a, b, c]))

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


def get_misc_notes(sv):
    """
    get frame and kinase domain information for
    fusion variants only
    sv -> str
    """
    misc_note = ""
    if "in frame" in sv.desc:
        frame_note = "is predicted to be in frame"
    else:
        frame_note = None
    kinase1, kinase2 = \
    map(lambda x: get_kinase_annotation(x, kinase_annotation), \
        (sv.fusionPartner1, sv.fusionPartner2))
    if kinase1 and kinase2:
        map(lambda x: kinase2.replace(x, ""),
            ("does not include the ", "includes "))
    misc_note = "and ".join(filter(lambda x: isinstance(x, str),\ 
        [frame_note, kinase1, kinase2]))
    if misc_note:
        return " This fusion %s." % misc_note
    else:
        return ""


def functional_significance(sv):
    """
    Determine if functional significance note is neccessary
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


def special_cases(sv):
    if sv.isIntragenic:
        if sv.bkp1.gene == "EGFR":
            re


special_case_notes = {
    "KDD" : "Note: The EGFR rearrangement is a kinase domain duplication (KDD).",
    "vIII" : "Note: The EGFR rearrangement is a vIII variant.",
    "CTD" : "Note: The EGFR rearrangement is a CTD variant.",
    "ERG" : "The fusion involves ERG non-canonical transcript (NM_00449)."
}
