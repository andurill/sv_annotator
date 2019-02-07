#!usr/bin/env python2

import os
import sys
import re
import warnings
import pandas as pd
import numpy as np
from main.models import bkp


def get_bkp_info(bkp, refFlat_summary, shift):
    """
    Get exon and intron features for a breakpoint object.
    bkp -> None
    """
    try:
        bkp_dict = refFlat_summary[
            (refFlat_summary['Gene'].values == bkp.gene) &
            (refFlat_summary['Transcript'].values == bkp.transcript)]\
            .to_dict('records').pop()
    except IndexError:
        raise Exception(
            "Exon and Intron information not found for %s in refFlat_summary."
            % (bkp.gene))
    #bkp.firstexon = bkp_dict['first_exon']
    # bkp.firstexon is hard-coded to 1 to follow
    # existing practice on IMPACT
    bkp.firstexon = "1"
    bkp.lastexon = str(bkp_dict['last_exon'])
    if bkp_dict['Strand'] == "+":
        bkp.startpos, bkp.stoppos = \
            bkp_dict['pos1'], bkp_dict['pos2']
    else:
        bkp.startpos, bkp.stoppos = \
            bkp_dict['pos2'], bkp_dict['pos1']
    bkp.exon, bkp.intron, bkp.site = [""]*3
    if bkp.desc.startswith("Exon "):
        bkp.exon = bkp.desc.split(" ")[1]
        bkp.site = "exon " + str(bkp.exon)
    elif bkp.desc.startswith("Intron"):
        exon = bkp.desc.split(" ")[5]
        if "before" in bkp.desc:
            bkp.intron = str(int(exon) - 1)
            if shift == 2:
                bkp.exon = exon
            else:
                bkp.exon = str(int(exon) - 1)
        else:
            bkp.intron = exon
            if shift == 1:
                bkp.exon = exon
            else:
                bkp.exon = str(int(exon) + 1)
        bkp.site = "intron " + str(bkp.intron)
    elif bkp.transcript == "NM_004449":
        bkp.exon = "4"
        bkp.site = "intron 3"
        warnings.warn("Breakpoint attributes are estimated for %s:%s. Manual review required!!!" % (
            bkp.chrom, bkp.pos), Warning)


def get_exon_order(bkp, order):
    if any([(bkp.strand == "+" and order == 1),
            (bkp.strand == "-" and order == 2),
            (order == 4)]):
        ordert = (bkp.exon, bkp.lastexon)
    elif any([(bkp.strand == "-" and order == 1),
              (bkp.strand == "+" and order == 2),
              (order == 3)]):
        ordert = (bkp.firstexon, bkp.exon)
    else:
        raise Exception(
            "Something terribly went wrong with \
                determining exon order for notes!")
    if ordert[0] == ordert[1]:
        return "exon %s" % (ordert[0])
    else:
        return "exons %s - %s" % ordert


def get_exons_involved(sv, refFlat_summary):
    """
    Get exons involved in an sv object based on the variant type
    and the breakpoint sites
    sv -> None
    """
    sv.bkpsites = ""
    if sv.isFusion:
        get_bkp_info(sv.fusionPartner1, refFlat_summary, 1)
        get_bkp_info(sv.fusionPartner2, refFlat_summary, 2)
        sv.bkpsites = get_bkpsite_note(
            sv, sv.fusionPartner1, sv.fusionPartner2)
        note1 = sv.fusionPartner1.gene + " " + \
            get_exon_order(sv.fusionPartner1, 3)
        note2 = sv.fusionPartner2.gene + " " + \
            get_exon_order(sv.fusionPartner2, 4)
        sv.fusionPartner1.variantSite1, sv.fusionPartner1.variantSite2 = \
            sv.fusionPartner1.startpos, sv.fusionPartner1.pos
        sv.fusionPartner2.variantSite1, sv.fusionPartner2.variantSite2 = \
            sv.fusionPartner2.pos, sv.fusionPartner2.stoppos
        if sv.isKnownFusion:
            return "%s and %s." % (note1, note2)
        else:
            return "of %s to %s." % (note1, note2)
    elif sv.bkp1.isPanel and sv.bkp2.isPanel and \
            sv.bkp1.isCoding and sv.bkp2.isCoding:
        get_bkp_info(sv.bkp1, refFlat_summary, 1)
        get_bkp_info(sv.bkp2, refFlat_summary, 2)
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s %s and %s %s" % \
                (sv.bkp1.gene, sv.bkp1.site,
                 sv.bkp2.gene, sv.bkp2.site)
            return "with breakpoints in %s." % (note1)
        elif sv.isIntragenic:
            get_bkp_info(sv.annotationPartner1, refFlat_summary, 2)
            get_bkp_info(sv.annotationPartner2, refFlat_summary, 1)
            sv.bkpsites = get_bkpsite_note(sv, sv.annotationPartner1,
                                           sv.annotationPartner2)
            if not sv.bkp1.exon == sv.bkp2.exon:
                note1 = "exons %s - %s" % (
                    sv.annotationPartner1.exon,
                    sv.annotationPartner2.exon)
            else:
                note1 = "exon %s" % (sv.annotationPartner1.exon)
            return "of %s." % (note1)
        else:
            sv.bkpsites = get_bkpsite_note(sv, sv.bkp1, sv.bkp2)
            note1 = sv.bkp1.gene + " " + get_exon_order(sv.bkp1, 1)
            note2 = sv.bkp2.gene + " " + get_exon_order(sv.bkp2, 2)
            return "of %s and %s." % (note1, note2)
    elif sv.bkp1.isPanel and sv.bkp1.isCoding:
        get_bkp_info(sv.bkp1, refFlat_summary, 1)
        if sv.svtype == "TRANSLOCATION":
            note1 = "with a breakpoint in %s" % (sv.bkp1.site)
        else:
            sv.bkpsites = get_bkpsite_note(sv, sv.bkp1, None)
            note1 = "of " + get_exon_order(sv.bkp1, 1)
        return "%s." % (note1)
    else:
        get_bkp_info(sv.bkp2, refFlat_summary, 2)
        if sv.svtype == "TRANSLOCATION":
            note1 = "with a breakpoint in %s" % (sv.bkp2.site)
        else:
            sv.bkpsites = get_bkpsite_note(sv, None, sv.bkp2)
            note1 = "of " + get_exon_order(sv.bkp2, 2)
        return "%s." % (note1)
    return


def getOverlap(a, b):
    # return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    return min(a[1], b[1]) - max(a[0], b[0])


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
    if sv.svtype == "INVERSION" or \
            sv.isIntragenic:
        conj = " is an "
    else:
        conj = " is a "
    if sv.isKnownFusion:
        prefix += str(sv.annotation.split(":")[0]) + \
            " involves "
    elif sv.isFusion:
        prefix += str(sv.annotation.split(":")[0]) + \
            conj + sv.svtype.lower() + \
            " that results in a fusion "
    elif sv.isIntragenic:
        prefix += str(sv.annotation.split(":")[0]) + \
            conj + "intragenic " + sv.svtype.lower() + " "
    else:
        prefix += str(sv.annotation.split(":")[0]) + \
            conj + sv.svtype.lower() + " "
    prefix = re.sub(r' \(NM_[0-9]+\)', "", prefix, count=2)
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
        return " Functional significance is undetermined."


def special_cases(sv):
    special_case_notes = {
        'KDD': 'Note: The EGFR rearrangement is a kinase domain duplication (KDD) alteration.',
        'vIII': 'Note: The EGFR rearrangement is a vIII alteration.',
        'CTD': 'Note: The EGFR rearrangement is a C-terminal domain (CTD) alteration.',
        'ERG': ' The structural variant involves the ERG non-canonical transcript (NM_004449).'
    }
    if sv.annotation.startswith("EGFR (NM_005228) rearrangement:") and \
            re.search(r"exons \d+ - \d+", sv.Note):
        se = re.search(r"exons \d+ - \d+", sv.Note)
        exon1, exon2 = se.group().replace("exons ","").split(" - ")
        if exon1 == "2" and exon2 == "7":
            return special_case_notes['vIII']
        elif exon1 in ("25", "26", "27", "28") and \
                exon2 in ("25", "26", "27", "28") and \
                sv.svtype == "DELETION":
            return special_case_notes['CTD']
        elif exon1 == "18" and \
            exon2 in ("25", "26") and \
                sv.svtype == "DUPLICATION":
            return special_case_notes['KDD']
        else:
            return sv.Note
    elif (sv.annotationPartner1.transcript == "NM_004449" and
          sv.annotationPartner1.isCoding) or \
            (sv.annotationPartner2.transcript == "NM_004449" and
                sv.annotationPartner2.isCoding):
        return sv.Note + special_case_notes['ERG']
    elif sv.annotationPartner1.gene == "CDKN2A" or \
            sv.annotationPartner2.gene == "CDKN2A":
        warnings.warn(
            "Manual review required for variants involving CDKN2A!!!", Warning)
        cdkn2a_tx = filter(lambda x: len(x) > 0,
                           set([sv.annotationPartner1.transcript,
                                sv.annotationPartner2.transcript]))
        if cdkn2a_tx.__len__() == 0:
            raise Exception("Cannot resolve CDKN2A isoforms.")
        elif cdkn2a_tx.__len__() == 2:
            cdkn2a_note = " This variant affects both CDKN2Ap14ARF (NM_058195) and CDKN2Ap16INK4A (NM_000077) isoforms."
        else:
            tx = cdkn2a_tx.pop()
            if tx == "NM_058195":
                cdkn2a_note = " This variant affects CDKN2Ap14ARF (NM_058195) isoform and may also affect CDKN2Ap16INK4A (NM_000077) isoform."
            else:
                cdkn2a_note = " This variant affects CDKN2Ap16INK4A (NM_000077) isoform and may also affect CDKN2Ap14ARF (NM_058195) isoform."
        return sv.Note + cdkn2a_note
    else:
        return sv.Note


def get_notes(sv, refFlat_summary, kinase_annotation):
    sv.prefix = get_prefix(sv)
    sv.exons = get_exons_involved(sv, refFlat_summary)
    if sv.isFusion:
        sv.misc = get_misc_notes(sv, kinase_annotation)
    else:
        sv.misc = ""
    sv.sig = (functional_significance(sv) or "")
    sv.Note = "".join([sv.prefix, sv.exons,
                       sv.bkpsites, sv.misc, sv.sig])
    position = get_position(sv)
    Note = special_cases(sv)
    return Note, position


def get_bkpsite_note(sv, b1=None, b2=None):
    if isinstance(b1, bkp) and b1.site.startswith("exon") and \
            isinstance(b2, bkp) and b2.site.startswith("exon"):
        if b1.gene == b2.gene:
            if b1.site == b2.site:
                bkps_note = " The breakpoints are within %s." % (b1.site)
            else:
                bkps_note = " The breakpoints are within %s and %s." % (
                    b1.site, b2.site)
        else:
            bkps_note = " The breakpoints are within %s %s and %s %s." % (
                b1.gene, b1.site, b2.gene, b2.site)
    else:
        bkps_note = " One of the breakpoints is within "
        if isinstance(b1, bkp) and b1.site.startswith("exon"):
            bkps_note += "%s %s." % (b1.gene,
                                     b1.site) if sv.isFusion else "%s." % (b1.site)
        elif isinstance(b2, bkp) and b2.site.startswith("exon"):
            bkps_note += "%s %s." % (b2.gene,
                                     b2.site) if sv.isFusion else "%s." % (b2.site)
        else:
            bkps_note = ""
    return bkps_note


def get_position(sv):
    note_local = sv.Note.split(".")[0]
    if sv.isFusion:
        position = "%s exon %s to %s exon %s" % \
            (sv.fusionPartner1.gene, sv.fusionPartner1.exon,
             sv.fusionPartner2.gene, sv.fusionPartner2.exon)
    else:
        position = re.sub(
            str(sv.prefix), "", note_local, count=1)
        position = re.sub(
            r'\bof \b|\binvolves \b|\bwith breakpoints in \b|\bwith a breakpoint in \b', "",
            position, count=1)
    return position
