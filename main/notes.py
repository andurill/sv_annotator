#!usr/bin/env python2

import os
import sys
import re
import pandas as pd
import numpy as np


def get_bkp_info(bkp, refFlat_summary):
    """
    Get exon and intron features for a breakpoint object.
    bkp -> None
    """
    bkp_dict = refFlat_summary[
        (refFlat_summary['Gene'].values == bkp.gene) &
        (refFlat_summary['Transcript'].values == bkp.transcript)]\
        .to_dict('records').pop()
    bkp.firstexon = bkp_dict['first_exon']
    bkp.lastexon = bkp_dict['last_exon']
    if bkp_dict['Strand'] == "+":
        bkp.startpos, bkp.stoppos = \
            bkp_dict['pos1'], bkp_dict['pos2']
    else:
        bkp.startpos, bkp.stoppos = \
            bkp_dict['pos2'], bkp_dict['pos1']
    bkp.exon, bkp.intron = [None]*2
    if bkp.desc.startswith("Exon "):
        bkp.exon = int(bkp.desc.split(" ")[1])
    elif bkp.desc.startswith("Intron"):
        exon = bkp.desc.split(" ")[5]
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


def get_exons_involved(sv, refFlat_summary):
    """
    Get exons involved in an sv object based on the variant type
    and the breakpoint sites
    sv -> None
    """
    if sv.isFusion:
        get_bkp_info(sv.fusionPartner1, refFlat_summary)
        get_bkp_info(sv.fusionPartner2, refFlat_summary)
        note1 = "%s exons %s - %s" % (sv.fusionPartner1.gene,
                                      sv.fusionPartner1.firstexon,
                                      sv.fusionPartner1.exon)
        note2 = "%s exons %s - %s" % (sv.fusionPartner2.gene,
                                      sv.fusionPartner2.exon,
                                      sv.fusionPartner2.lastexon)
        #sv.fusionPartner1.variantSite1 = sv.fusionPartner1.firstexon
        #sv.fusionPartner1.variantSite2 = sv.fusionPartner1.exon
        #sv.fusionPartner2.variantSite1 = sv.fusionPartner2.firstexon
        #sv.fusionPartner2.variantSite2 = sv.fusionPartner2.exon
        sv.fusionPartner1.variantSite1, sv.fusionPartner1.variantSite2 = \
            sv.fusionPartner1.startpos, sv.fusionPartner1.pos
        sv.fusionPartner2.variantSite1, sv.fusionPartner2.variantSite2 = \
            sv.fusionPartner2.pos, sv.fusionPartner2.stoppos
        if sv.isKnownFusion:
            return "%s and %s." % (note1, note2)
        else:
            return "%s to %s." % (note1, note2)
    elif sv.bkp1.isPanel and sv.bkp2.isPanel and \
            sv.bkp1.isCoding and sv.bkp2.isCoding:
        get_bkp_info(sv.annotationPartner1, refFlat_summary)
        get_bkp_info(sv.annotationPartner2, refFlat_summary)
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s %s and %s %s" % (sv.annotationPartner1.gene,
                                         sv.annotationPartner1.site,
                                         sv.annotationPartner2.gene,
                                         sv.annotationPartner2.site)
            return "%s." % (note1)
        elif sv.isIntragenic:
            note1 = "exons %s - %s" % (sv.annotationPartner1.exon,
                                       sv.annotationPartner2.exon)
            return "%s." % (note1)
        else:
            note1 = "%s exons %s - %s" % (sv.annotationPartner1.gene,
                                          sv.annotationPartner1.firstexon,
                                          sv.annotationPartner1.exon)
            note2 = "%s exons %s - %s" % (sv.annotationPartner2.gene,
                                          sv.annotationPartner2.exon,
                                          sv.annotationPartner2.lastexon)
            return "%s and %s." % (note1, note2)
    elif sv.bkp1.isPanel and sv.bkp1.isCoding:
        get_bkp_info(sv.annotationPartner1, refFlat_summary)
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s" % (sv.annotationPartner1.site)
        elif sv.bkp1.strand == "+":
            note1 = "exons %s - %s" % (sv.annotationPartner1.exon,
                                       sv.annotationPartner1.lastexon)
        elif sv.bkp1.strand == "-":
            note1 = "exons %s - %s" % (sv.annotationPartner1.firstexon,
                                       sv.annotationPartner1.exon)
        return "%s." % (note1)
    else:
        get_bkp_info(sv.annotationPartner2, refFlat_summary)
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s" % (sv.annotationPartner2.site)
        elif sv.bkp2.strand == "-":
            note1 = "exons %s - %s" % (sv.annotationPartner2.exon,
                                       sv.annotationPartner2.lastexon)
        elif sv.bkp2.strand == "+":
            note1 = "exons %s - %s" % (sv.annotationPartner2.firstexon,
                                       sv.annotationPartner2.exon)
        return "%s." % (note1)
    return


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_kinase_annotation(bkp, kinase_annotation):
    bkp.isEntireKinase = False
    kinase_interval = []
    try:
        interval = [bkp.variantSite1, bkp.variantSite2]
        if bkp.strand == '-':
            interval.sort()
        kinase_interval = kinase_annotation[
            (kinase_annotation['HUGO'].values == bkp.gene)].\
            iloc[0, [1, 4]].values.tolist()
    except IndexError:
        pass
        # no values found
    except NameError:
        pass
        # not fusion
    except ValueError:
        pass
        # what
    if kinase_interval:
        interval, kinase_interval = \
            map(lambda x: [int(y) for y in x], (interval, kinase_interval))
        overlap = getOverlap(interval, kinase_interval)
        if overlap == 0:
            return "does not include the kinase domain of %s" % (bkp.gene)
        else:
            if min(interval) <= min(kinase_interval) and \
                    max(interval) >= max(kinase_interval):
                bkp.isEntireKinase = True
                return "includes the kinase domain of %s" % (bkp.gene)
            else:
                return "includes a part of the kinase domain of %s" % (bkp.gene)
    return None


def get_prefix(sv):
    """
    Get the prefix of a note based on variant type
    sv -> str
    """
    prefix = "The "
    fusion_type = ""  # "rearrangement"
    if sv.isKnownFusion:
        fusion_type = ""  # "fusion"
        prefix += str(sv.annotation.split(":")[0]) + \
            fusion_type + " involves "
    elif sv.isFusion:
        prefix += str(sv.annotation.split(":")[0]) + \
            fusion_type + " is a " + sv.svtype.lower() + \
            " that results in a fusion of "
    elif sv.isIntragenic:
        prefix += str(sv.annotation.split(":")[0]) + \
            fusion_type + " is an intragenic " + \
            sv.svtype.lower() + " of "
    else:
        prefix += str(sv.annotation.split(":")[0]) + \
            fusion_type + " is a " + sv.svtype.lower() + " of "
    return prefix


def get_misc_notes(sv, kinase_annotation):
    """
    get frame and kinase domain information for
    fusion variants only
    sv -> str
    """
    misc_note = ""
    if "in frame" in sv.description:
        frame_note = "is predicted to be in frame"
    else:
        frame_note = None
    kinase1, kinase2 = \
        map(lambda x: get_kinase_annotation(x, kinase_annotation),
            (sv.fusionPartner1, sv.fusionPartner2))
    if kinase1 and kinase2:
        map(lambda x: kinase2.replace(x, ""),
            ("does not include the ", "includes "))
    misc_note = " and ".join(filter(lambda x: isinstance(x, str),
                                    [frame_note, kinase1, kinase2]))
    if misc_note:
        return " The fusion %s." % (misc_note)
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


def special_cases(sv, Note):
    special_case_notes = {
        'KDD': 'Note: The EGFR rearrangement is a kinase domain duplication (KDD).',
        'vIII': 'Note: The EGFR rearrangement is a vIII variant.',
        'CTD': 'Note: The EGFR rearrangement is a CTD variant.',
        'ERG': ' The structural variant involves the ERG non-canonical transcript (NM_004449).'
    }
    if sv.isIntragenic:
        if sv.annotationPartner1.gene == "EGFR":
            if sv.annotationPartner1.exon == 2 and \
                    sv.annotationPartner2.exon == 7:
                return special_case_notes['vIII']
        else:
            return Note
    elif (sv.bkp1.transcript == "NM_004449" and sv.bkp1.isCoding) or \
            (sv.bkp2.transcript == "NM_004449" and sv.bkp2.isCoding):
        return Note + special_case_notes['ERG']
    else:
        return Note


def get_notes(sv, refFlat_summary, kinase_annotation):
    Note = re.sub(r' \(NM_[0-9]+\)', "", get_prefix(sv), count=2)
    Note += get_exons_involved(sv, refFlat_summary)
    if sv.isFusion:
        Note += get_misc_notes(sv, kinase_annotation)
    Note += functional_significance(sv)
    position = get_position(Note)
    Note = special_cases(sv, Note)
    return Note


def get_position(sv, Note):
    if sv.isFusion:
        return "%s %s to %s %s" % (sv.fusionPartner1.gene,
                                   sv.fusionPartner1.exon, 
                                   sv.fusionPartner2.gene,
                                   sv.fusionPartner2.exon)
    else:
        position = re
