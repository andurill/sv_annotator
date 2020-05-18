#!/usr/bin/env python2

import os
import sys
import re
import logging
import pandas as pd
import numpy as np
from main.models import bkp

logger = logging.getLogger("basic_logger")


def get_bkp_info(bkp, refFlat_summary, orientation, fusion=0):
    """
    Get exon and intron features for a breakpoint object
    bkp, df, int, int -> None
    """
    try:
        bkp_dict = (
            refFlat_summary[
                (refFlat_summary["Gene"].values == bkp.gene)
                & (refFlat_summary["Transcript"].values == bkp.transcript)
            ]
            .to_dict("records")
            .pop()
        )
    except IndexError:
        raise Exception(
            "Exon and Intron information not found for %s in refFlat_summary."
            % (bkp.gene)
        )
    # bkp.firstexon = bkp_dict['first_exon']
    # bkp.firstexon is hard-coded to 1 to follow
    # existing practice on IMPACT
    bkp.firstexon = "1"
    bkp.lastexon = str(bkp_dict["last_exon"])
    if bkp_dict["Strand"] == "+":
        bkp.startpos, bkp.stoppos = bkp_dict["pos1"], bkp_dict["pos2"]
    else:
        bkp.startpos, bkp.stoppos = bkp_dict["pos2"], bkp_dict["pos1"]
    bkp.exon, bkp.intron, bkp.site = [""] * 3
    if bkp.desc.startswith("Exon "):
        bkp.exon = bkp.desc.split(" ")[1]
        bkp.site = "exon " + str(bkp.exon)
        if get_bkp_type(bkp, fusion, orientation) == 1:
            bkp.variantSite1, bkp.variantSite2 = bkp.pos, bkp.stoppos
        else:
            bkp.variantSite1, bkp.variantSite2 = bkp.startpos, bkp.pos
    elif bkp.desc.startswith("Intron "):
        exon = bkp.desc.split(" ")[-1]
        if "after" in bkp.desc:
            bkp.intron = exon
            if get_bkp_type(bkp, fusion, orientation) == 1:
                bkp.exon, bkp.variantSite1, bkp.variantSite2 = (
                    str(int(exon) + 1),
                    bkp.pos,
                    bkp.stoppos,
                )
            else:
                bkp.exon, bkp.variantSite1, bkp.variantSite2 = (
                    exon,
                    bkp.startpos,
                    bkp.pos,
                )
        elif "before" in bkp.desc:
            bkp.intron = str(int(exon) - 1)
            if get_bkp_type(bkp, fusion, orientation) == 2:
                bkp.exon, bkp.variantSite1, bkp.variantSite2 = (
                    str(int(exon) - 1),
                    bkp.startpos,
                    bkp.pos,
                )
            else:
                bkp.exon, bkp.variantSite1, bkp.variantSite2 = (
                    exon,
                    bkp.pos,
                    bkp.stoppos,
                )
        bkp.site = "intron " + bkp.intron
    elif bkp.transcript == "NM_004449" and fusion == 1:
        bkp.exon, bkp.site, bkp.variantSite1, bkp.variantSite2 = (
            "4",
            "intron 3",
            bkp.pos,
            bkp.stoppos,
        )
        logger.warning(
            "Breakpoint attributes are estimated for %s:%s. Manual review required"
            % (bkp.chrom, bkp.pos)
        )


def get_bkp_type(bkp, fusion, orientation):
    if any(
        [
            fusion == 1 and orientation == 2,
            fusion == 0 and bkp.strand == "+" and orientation == 1,
            fusion == 0 and bkp.strand == "-" and orientation == 2,
        ]
    ):
        return 1
    elif any(
        [
            fusion == 1 and orientation == 1,
            fusion == 0 and bkp.strand == "+" and orientation == 2,
            fusion == 0 and bkp.strand == "-" and orientation == 1,
        ]
    ):
        return 2
    else:
        raise Exception("Something went terribly wrong with breakpoint expansion.")


def get_exon_order(bkp, order):
    """
    Determine gene order for notes
    bkp, int -> str
    """
    if any(
        [
            (bkp.strand == "+" and order == 1),
            (bkp.strand == "-" and order == 2),
            (order == 4),
        ]
    ):
        ordert = (bkp.exon, bkp.lastexon)
    elif any(
        [
            (bkp.strand == "-" and order == 1),
            (bkp.strand == "+" and order == 2),
            (order == 3),
        ]
    ):
        ordert = (bkp.firstexon, bkp.exon)
    else:
        raise Exception(
            "Something terribly went wrong with \
                determining exon order for notes!"
        )
    if ordert[0] == ordert[1]:
        return "exon %s" % (ordert[0])
    else:
        return "exons %s - %s" % ordert


def get_exons_involved(sv, refFlat_summary):
    """
    Get exons involved in an sv object based on the variant type
    and the breakpoint sites
    (sv, df) -> None
    """
    sv.bkpsites = ""
    if sv.isFusion:
        get_bkp_info(sv.fusionPartner1, refFlat_summary, 1, 1)
        get_bkp_info(sv.fusionPartner2, refFlat_summary, 2, 1)
        sv.bkpsites = get_bkpsite_note(sv, sv.fusionPartner1, sv.fusionPartner2)
        note1 = sv.fusionPartner1.gene + " " + get_exon_order(sv.fusionPartner1, 3)
        note2 = sv.fusionPartner2.gene + " " + get_exon_order(sv.fusionPartner2, 4)
        if sv.isKnownFusion:
            sv.exons = "%s and %s." % (note1, note2)
        else:
            sv.exons = "%s to %s." % (note1, note2)
    elif all(
        [
            sv.annotationPartner1.isPanel,
            sv.annotationPartner2.isPanel,
            sv.annotationPartner1.isCoding,
            sv.annotationPartner2.isCoding,
        ]
    ):
        get_bkp_info(sv.annotationPartner1, refFlat_summary, 1)
        get_bkp_info(sv.annotationPartner2, refFlat_summary, 2)
        if sv.svtype == "TRANSLOCATION":
            sv.exons = "%s %s and %s %s." % (
                sv.annotationPartner1.gene,
                sv.annotationPartner1.site,
                sv.annotationPartner2.gene,
                sv.annotationPartner2.site,
            )
        elif sv.isIntragenic:
            intra1, intra2 = (1, 2) if sv.annotationPartner1.strand == "+" else (2, 1)
            get_bkp_info(sv.annotationPartner1, refFlat_summary, intra1)
            get_bkp_info(sv.annotationPartner2, refFlat_summary, intra2)
            sv.bkpsites = get_bkpsite_note(
                sv, sv.annotationPartner1, sv.annotationPartner2
            )
            if all(
                [
                    sv.annotationPartner1.intron,
                    sv.annotationPartner2.intron,
                    sv.annotationPartner1.intron == sv.annotationPartner2.intron,
                ]
            ):
                sv.exons = "intron %s." % (sv.annotationPartner1.intron)
            elif not sv.annotationPartner1.exon == sv.annotationPartner2.exon:
                sv.exons = "exons %s - %s." % (
                    sv.annotationPartner1.exon,
                    sv.annotationPartner2.exon,
                )
            else:
                sv.exons = "exon %s." % (sv.annotationPartner1.exon)
            sv.annotationPartner1.variantSite1, sv.annotationPartner2.variantSite1 = [
                sv.annotationPartner1.pos
            ] * 2
            sv.annotationPartner1.variantSite2, sv.annotationPartner2.variantSite2 = [
                sv.annotationPartner2.pos
            ] * 2
        else:
            sv.bkpsites = get_bkpsite_note(
                sv, sv.annotationPartner1, sv.annotationPartner2
            )
            note1 = (
                sv.annotationPartner1.gene
                + " "
                + get_exon_order(sv.annotationPartner1, 1)
            )
            note2 = (
                sv.annotationPartner2.gene
                + " "
                + get_exon_order(sv.annotationPartner2, 2)
            )
            sv.exons = "%s and %s." % (note1, note2)
    elif sv.annotationPartner1.isPanel and sv.annotationPartner1.isCoding:
        get_bkp_info(sv.annotationPartner1, refFlat_summary, 1)
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s" % (sv.annotationPartner1.site)
        else:
            sv.bkpsites = get_bkpsite_note(sv, sv.annotationPartner1, None)
            note1 = get_exon_order(sv.annotationPartner1, 1)
        sv.exons = "%s." % (note1)
    else:
        get_bkp_info(sv.annotationPartner2, refFlat_summary, 2)
        if sv.svtype == "TRANSLOCATION":
            note1 = "%s" % (sv.annotationPartner2.site)
        else:
            sv.bkpsites = get_bkpsite_note(sv, None, sv.annotationPartner2)
            note1 = get_exon_order(sv.annotationPartner2, 2)
        sv.exons = "%s." % (note1)
    return


def getOverlap(a, b):
    """
    Determine if two intervals overlap
    list, list -> int
    """
    a.sort(), b.sort()
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_kinase_status(bkp, kinase_annotation):
    """
    Get kinase domain annotation for a given bkp
    bkp, df -> (str or None)
    """
    bkp.isEntireKinase = None
    if not bkp.isKinase:
        return
    interval, kinase_interval = [[]] * 2
    try:
        interval.extend([bkp.variantSite1, bkp.variantSite2])
        if bkp.strand == "-":
            interval.sort()
        kinase_interval = (
            kinase_annotation[(kinase_annotation["HUGO"].values == bkp.gene)]
            .iloc[0, [1, 4]]
            .values.tolist()
        )
        kinase_interval.sort()
    except AttributeError:
        e = (
            "Breakpoint site positions cannot be determined "
            "for gene %s. Check if the raw SV call was made on "
            "Canonical transcript."
        ) % (bkp.gene)
        logger.warning(e)
    except (IndexError, NameError, ValueError):
        logger.warning("Cannot find kinase domain information for %s" % (bkp.gene))
    if interval and kinase_interval:
        interval, kinase_interval = map(
            lambda x: [int(y) for y in x], (interval, kinase_interval)
        )
        overlap = getOverlap(interval, kinase_interval)
        if overlap == 0:
            bkp.isEntireKinase = 0
        else:
            if min(interval) <= min(kinase_interval) and max(interval) >= max(
                kinase_interval
            ):
                bkp.isEntireKinase = 1
            else:
                bkp.isEntireKinase = -1
    return


def get_kinase_note(bkp):
    """
    Get kinase notes based on kinase status of the bkp
    bkp -> (str or None)
    """
    if not bkp.isKinase:
        return None
    else:
        if bkp.isEntireKinase == 1:
            return "includes the kinase domain of %s" % (bkp.gene)
        elif bkp.isEntireKinase == -1:
            return "includes a part of the kinase domain of %s" % (bkp.gene)
        else:
            return "does not include the kinase domain of % s" % (bkp.gene)


def get_prefix(sv):
    """
    Get the prefix of a note based on variant type
    sv -> None
    """
    prefix = "The "
    if sv.svtype == "INVERSION" or sv.isIntragenic:
        conj = " is an "
    else:
        conj = " is a "
    if sv.isKnownFusion:
        prefix += str(sv.annotation.split(":")[0]) + " involves "
    elif sv.isFusion:
        prefix += (
            str(sv.annotation.split(":")[0])
            + conj
            + sv.svtype.lower()
            + " that results in a fusion of "
        )
    elif sv.isIntragenic:
        if all(
            [
                sv.annotationPartner1.intron,
                sv.annotationPartner2.intron,
                sv.annotationPartner1.intron == sv.annotationPartner2.intron,
            ]
        ):
            prefix += (
                str(sv.annotation.split(":")[0])
                + conj
                + "intragenic "
                + sv.svtype.lower()
                + " with breakpoints in "
            )
        else:
            prefix += (
                str(sv.annotation.split(":")[0])
                + conj
                + "intragenic "
                + sv.svtype.lower()
                + " of "
            )
    elif sv.svtype == "TRANSLOCATION":
        if all(
            [
                sv.annotationPartner1.isPanel,
                sv.annotationPartner2.isPanel,
                sv.annotationPartner1.isCoding,
                sv.annotationPartner2.isCoding,
            ]
        ):
            TRA_join = " with breakpoints in "
        else:
            TRA_join = " with a breakpoint in "
        prefix += str(sv.annotation.split(":")[0]) + conj + sv.svtype.lower() + TRA_join
    else:
        prefix += str(sv.annotation.split(":")[0]) + conj + sv.svtype.lower() + " of "
    sv.prefix = re.sub(r" \(NM_[0-9]+\)", "", prefix, count=2)


def get_misc_notes(sv):
    """
    get frame and kinase domain information for
    fusion variants only
    sv -> None
    """
    if sv.svtype == "TRANSLOCATION" and not sv.isFusion:
        sv.misc = ""
        return
    prefix = "The fusion" if sv.isFusion else "The rearrangement"
    misc_note, frame_note = [None] * 2
    if sv.isFusion and "in frame" in sv.description:
        frame_note = "is predicted to be in frame"
    if sv.isIntragenic:
        kinase1, kinase2 = get_kinase_note(sv.annotationPartner1), None
    else:
        kinase1, kinase2 = map(
            lambda x: get_kinase_note(x), (sv.annotationPartner1, sv.annotationPartner2)
        )
    if all(
        [
            kinase1,
            kinase2,
            sv.annotationPartner1.isEntireKinase
            == sv.annotationPartner2.isEntireKinase,
        ]
    ):
        regex_str = (
            r"\b"
            + "does not include the kinase domain of "
            + r"\b|"
            + r"\b"
            + "includes the kinase domain of "
            + r"\b|"
            + r"\b"
            + "includes a part of the kinase domain of "
            + r"\b"
        )
        kinase2 = re.sub(regex_str, "", kinase2, count=1)
        kinase1 = kinase1.replace("domain", "domains").replace("a part", "parts")
    misc_note = " and ".join(
        filter(lambda x: isinstance(x, str), [frame_note, kinase1, kinase2])
    )
    if misc_note:
        sv.misc = " %s %s." % (prefix, misc_note)
    else:
        sv.misc = ""


def clinical_warning_note(sv):
    bkps_involved = filter(
        lambda x: x.isKinase and x.isHotspot,
        [sv.annotationPartner1, sv.annotationPartner2],
    )
    bkps_genes = list(set(map(lambda x: x.gene, bkps_involved)))
    genes_involved = " and ".join(bkps_genes)
    warning_note = (
        " Functional significance is undetermined and "
        "further testing to determine the presence "
        "or absence of a targetable oncogenic fusion involving %s "
        "is required. This sample has been nominated for further "
        "analysis using the Archer targeted RNAseq assay. "
        "Archer will be performed and reported under a separate "
        "accession number if additional material is available."
    ) % (genes_involved)
    # basic_warning_note not used in the current version
    basic_warning_note = (
        " Functional significance is undetermined. "
        "Further testing to characterize the structural variant "
        "involving %s is required."
    ) % (genes_involved)
    return warning_note


def functional_significance(sv):
    """
    Determine if functional significance note is neccessary
    based on variant type
    sv -> None
    """
    func_sig_note = ""
    if sv.isKnownFusion:
        pass
    elif all([sv.annotationPartner1.isKinase, sv.annotationPartner1.isHotspot]) or all(
        [sv.annotationPartner2.isKinase, sv.annotationPartner2.isHotspot]
    ):
        func_sig_note = clinical_warning_note(sv)
    elif sv.svtype == "DELETION" and any(
        [
            all([sv.bkp1.isPanel, sv.bkp1.isCoding, sv.bkp1.isTumourSuppressor]),
            all([sv.bkp2.isPanel, sv.bkp2.isCoding, sv.bkp2.isTumourSuppressor]),
        ]
    ):
        pass
    else:
        func_sig_note = " Functional significance is undetermined."
    sv.sig = func_sig_note


def special_cases(sv):
    """
    Returns a customized note for special cases
    sv -> str
    """
    special_case_notes = {
        "KDD": "The EGFR rearrangement is a kinase domain duplication (KDD) alteration.",
        "vIII": "The EGFR rearrangement is a vIII alteration.",
        "CTD": "The EGFR rearrangement is a C-terminal domain (CTD) alteration.",
        "ERG": "The structural variant involves the ERG non-canonical transcript (NM_004449).",
    }
    custom_note = sv.Note
    if sv.annotation.startswith("EGFR (NM_005228) rearrangement:"):
        if re.search(r"exons \d+ - \d+", sv.Note):
            exon1, exon2 = (
                re.search(r"exons \d+ - \d+", sv.Note)
                .group()
                .replace("exons ", "")
                .split(" - ")
            )
            if all([exon1 == "2", exon2 == "7", sv.svtype != "DUPLICATION"]):
                custom_note = special_case_notes["vIII"]
            elif all(
                [
                    exon1 in ("25"),
                    exon2 in ("26", "27", "28"),
                    sv.svtype != "DUPLICATION",
                ]
            ):
                custom_note = special_case_notes["CTD"]
            elif all(
                [exon1 == "18", exon2 in ("25", "26"), sv.svtype == "DUPLICATION"]
            ):
                custom_note = special_case_notes["KDD"]
    elif (
        sv.annotationPartner1.transcript == "NM_004449"
        and sv.annotationPartner1.isCoding
    ) or (
        sv.annotationPartner2.transcript == "NM_004449"
        and sv.annotationPartner2.isCoding
    ):
        custom_note += " " + special_case_notes["ERG"]
    elif "CDKN2A" in (sv.annotationPartner1.gene, sv.annotationPartner2.gene):
        logger.warning("Manual review required for variants involving CDKN2A!!!")
        cdkn2a_tx = filter(
            lambda x: isinstance(x, str),
            set([sv.annotationPartner1.transcript, sv.annotationPartner2.transcript]),
        )
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
        custom_note += cdkn2a_note
    return custom_note


def get_bkpsite_note(sv, b1=None, b2=None):
    """
    Determine breakpoint details.
    sv, (None or bkp), (None or bkp) -> str
    """
    if (
        isinstance(b1, bkp)
        and b1.site.startswith("exon")
        and isinstance(b2, bkp)
        and b2.site.startswith("exon")
    ):
        if b1.gene == b2.gene:
            if b1.site == b2.site:
                bkps_note = " The breakpoints are within %s." % (b1.site)
            else:
                bkps_note = " The breakpoints are within %s and %s." % (
                    b1.site,
                    b2.site,
                )
        else:
            bkps_note = " The breakpoints are within %s %s and %s %s." % (
                b1.gene,
                b1.site,
                b2.gene,
                b2.site,
            )
    else:
        bkps_note = " One of the breakpoints is within "
        if isinstance(b1, bkp) and b1.site.startswith("exon"):
            bkps_note += (
                "%s %s." % (b1.gene, b1.site) if sv.isFusion else "%s." % (b1.site)
            )
        elif isinstance(b2, bkp) and b2.site.startswith("exon"):
            bkps_note += (
                "%s %s." % (b2.gene, b2.site) if sv.isFusion else "%s." % (b2.site)
            )
        else:
            bkps_note = ""
    return bkps_note


def get_position(sv):
    """
    Derive position string from notes
    sv -> str
    """
    note_local = sv.Note.split(".")[0]
    if sv.isFusion:
        position = "%s exon %s to %s exon %s" % (
            sv.fusionPartner1.gene,
            sv.fusionPartner1.exon,
            sv.fusionPartner2.gene,
            sv.fusionPartner2.exon,
        )
    else:
        position = re.sub(str(sv.prefix), "", note_local, count=1)
        position = re.sub(
            r"\bof \b|\binvolves \b|\bwith breakpoints in \b|\bwith a breakpoint in \b",
            "",
            position,
            count=1,
        )
    return position


def override_fusion(sv):
    if not sv.isFusion:
        pass
    else:
        fusion_type1, fusion_type2 = "rearrangement", "fusion"
        fusion_conj1, fusion_conj2 = "to", "and"
        if any(
            [
                all(
                    [
                        sv.fusionPartner1.isKinase,
                        sv.fusionPartner1.isHotspot,
                        sv.fusionPartner1.isEntireKinase == 1,
                    ]
                ),
                all(
                    [
                        sv.fusionPartner2.isKinase,
                        sv.fusionPartner2.isHotspot,
                        sv.fusionPartner2.isEntireKinase == 1,
                    ]
                ),
            ]
        ):
            sv.isKnownFusion = True
            sv.annotation = sv.annotation.replace(fusion_type1, fusion_type2)
            sv.exons = sv.exons.replace(fusion_conj1, fusion_conj2)

        elif any(
            [
                all(
                    [
                        sv.fusionPartner1.isKinase,
                        sv.fusionPartner1.isHotspot,
                        sv.fusionPartner1.isEntireKinase != 1,
                    ]
                ),
                all(
                    [
                        sv.fusionPartner2.isKinase,
                        sv.fusionPartner2.isHotspot,
                        sv.fusionPartner2.isEntireKinase != 1,
                    ]
                ),
            ]
        ):
            sv.isKnownFusion = False
            sv.annotation = sv.annotation.replace(fusion_type2, fusion_type1)
            sv.exons = sv.exons.replace(fusion_conj2, fusion_conj1)
    return


def get_sv_oncokb_type(sv):
    """
    This function returns the oncokb sv type `FUSION`, `vIII`, `KDD` or `CTD` 
    sv -> str
    """
    if sv.isKnownFusion:
        return "FUSION"
    else:
        special_vars = ["vIII", "KDD", "CTD"]
        for i, var in enumerate(special_vars):
            if var in sv.Note and "EGFR" in sv.Note:
                return special_vars[i]


def get_notes(sv, refFlat_summary, kinase_annotation):
    """
    Main note function to call relevant helper functions
    sv, df, df -> tuple
    """
    # Get exons and breakpoints invovled in SV
    get_exons_involved(sv, refFlat_summary)
    # Get kinase domain annotation
    map(lambda bkp: get_kinase_status(bkp, kinase_annotation), (sv.bkp1, sv.bkp2))
    # Override status of known fusion if necessary
    override_fusion(sv)
    # Get the first statement of clinical SV note
    get_prefix(sv)
    # Get miscellaneous notes
    get_misc_notes(sv)
    # Determine the clinical significance of the SV
    functional_significance(sv)
    sv.Note = "".join([sv.prefix, sv.exons, sv.bkpsites, sv.misc, sv.sig])
    position = get_position(sv)
    sv.Note = ": ".join(["Note", special_cases(sv)])
    oncokb_sv_type = get_sv_oncokb_type(sv)
    return sv.Note, position, oncokb_sv_type
