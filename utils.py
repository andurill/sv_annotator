'''
Copyright (c) 2019 Memorial Sloan-Kettering Cancer Center.

This script is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 2.1 of the License, or
any later version.

This script is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
documentation provided hereunder is on an "as is" basis, and
Memorial Sloan-Kettering Cancer Center has no obligations to provide
maintenance, support, updates, enhancements or modifications.  In no
event shall Memorial Sloan-Kettering Cancer Center be liable to any
party for direct, indirect, special, incidental or consequential damages,
including lost profits, arising out of the use of this software and its
documentation, even if Memorial Sloan-Kettering Cancer Center has been
advised of the possibility of such damage.  See the GNU Lesser General
Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; if not, write to the Free Software Foundation,
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.

Created on Jan 20, 2019

@author: balakra1@mskcc.org,jayakumg@mskcc.org
@version: 1.0.0

'''


import requests, re
from constants.constants import *

def make_get_request(query,logger):
    """
    Generate query string in the format required by vep
    str -> request
    """
    request = None
    __server_url = "{}{}".format(EMSAMBLE_SERVER,query)
    try:
        request = requests.get( __server_url, headers=API_SET_CONTENT_TYPE_JSON, timeout=TIMEOUT)
    except requests.exceptions.RequestException as e:
        logger.error("RequestException while querying vep for {}".format(query))
        logger.error(e)
    except Exception as e:
        logger.error("Exception while querying vep for {}".format(query))
        logger.error(e)
    finally:
        return request

def make_query(bkp, dummy_ref):
    """
    Make dummy query string for the second
    request sent to vep server
    bkp, str -> str
    """
    chrom, coord = bkp.get_chrom(), bkp.get_pos()
    revcomp = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N"
    }
    dummy_alt = revcomp[dummy_ref]
    query = "{chrom}:g.{coord}{dummy_ref}>{dummy_alt}?".format(
                                                                chrom=chrom,
                                                                coord=str(coord),
                                                                dummy_ref=dummy_ref,
                                                                dummy_alt=dummy_alt
                                                            )
    return query


def get_cdna_pos(bkp,logger):
    """
    Get cdna position for a bkp object by querying
    vep server
    bkp -> tuple
    """
    dummy_ref = "C"
    query = make_query(bkp, dummy_ref)
    cdna, tx = [""]*2
    try:
        request = make_get_request(query,logger)

        if not request.ok and not request:
            s = request.text
            actual_ref = s[s.find("(")+1:s.find(")")]
            if actual_ref not in ["A", "C", "G", "T"] or \
                    actual_ref == dummy_ref:
                request.raise_for_status()
            else:
                query = make_query(bkp, actual_ref)
                request = make_get_request(query,logger)
                if not request.ok and not request:
                    request.raise_for_status()

        decoded = request.json()
        try:
            result = dict(str(s).split(':', 1) for s in decoded[0]['hgvsc'])
        except KeyError:
            result = {}
            logger.error("Cannot find any cDNA annotations for gene {}".format(bkp.get_gene()))
     
        for tx in bkp.get_transcript():
            if tx in result:
                query_cdna = result[tx]
                if not (query_cdna.startswith("c.-") or query_cdna.startswith("c.*") ):
                    cdna = query_cdna[:-3]
                    break

        if cdna == "":
            cdna = "chr" + bkp.get_chrom() + ":g." + str(bkp.get_pos())
    except Exception as e:
        logger.error("Exception in get_cdna_pos")
        logger.error(e)
    return tx, cdna

def is_oncokb_fusion(__sv,logger):
    r"""
    THs fucntion queries the oncokb API and set the __is_known_oncokb_fusion flag

    Parameters
    ----------
    __sv : SVAnnotator
    logger : loggerAdapter

    Returns
    -------
    __is_known_oncokb_fusion : bool

    Notes
    ---------
    __is_known_oncokb_fusion is set to true only if getOncogenicStatus() from oncokb API is `oncogenic`

    """ 
    __bk1 = __sv.get_annotation_partner_1()
    __bk2 = __sv.get_annotation_partner_2()

    __fusion_genes = "{}:{}".format(__bk1.get_gene(),__bk2.get_gene())
    return True if __fusion_genes in ONCOKB_FUSIONS else False

def get_cytoband(__bkp,logger):
    """
    Get cytoband of a breakpoint by querying panda dataframe
    bkp -> str
    """
    __cytoband = CB_DF[(CB_DF['Chr'].values == __bkp.get_chrom()) &
                           (CB_DF['Bp_start'].values <= __bkp.get_pos()) &
                           (CB_DF['Bp_stop'].values >= __bkp.get_pos())]
    if len(__cytoband) == 0:
        logger.error("Cytobad Missing")
        raise ("Cytobad Missing")
    elif len(__cytoband) > 1:
        logger.error("Mutliple Cytobad")
        raise ("Mutliple Cytobad")
    locus = __cytoband.values[0]
    return locus[0] + locus[1]

def  get_coordinate(a1,a2,logger):
    __cband_1 = get_cytoband(a1,logger)
    __cband_2 = get_cytoband(a2,logger)
    __coordinate = "t({chr1};{chr2})({cb1};{cb2})(chr{chr1}:g.{pos1}::chr{chr2}:g.{pos2})".format(chr1=a1.get_chrom(),chr2=a2.get_chrom(),cb1=__cband_1,cb2=__cband_2,pos1=a1.get_pos(),pos2=a2.get_pos())
    return __coordinate

def get_overlap(a, b):
    """
    Determine if two intervals overlap
    list, list -> int
    """
    a.sort(), b.sort()
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def get_bkp_info(bkp, orientation=0, fusion=0, logger=None):
    """
    Get exon and intron features for a breakpoint object.
    bkp -> None
    """
    try:
        try:
            bkp_dict = REFFLAT_SUMMARY[(REFFLAT_SUMMARY['Gene'].values == bkp.get_gene()) & (REFFLAT_SUMMARY['Transcript'].values == bkp.get_transcript())].to_dict('records').pop()
        except IndexError:
            logger.error("Exon and Intron information not found for %s in refFlat_summary."% (bkp.get_gene()))
            raise Exception("Exon and Intron information not found for %s in refFlat_summary."% (bkp.get_gene()))

        bkp.set_last_exon(bkp_dict['last_exon'])
        bkp.set_start_pos(bkp_dict['pos2'])
        bkp.set_stop_pos(bkp_dict['pos1'])
        __bkp_type = get_bkp_type(bkp ,fusion, orientation)

        if __bkp_type == 0 : raise Exception("Something went wrong with bkp expamnsion")
        
        if bkp_dict['Strand'] == POSITIVE_STRAND:
            bkp.set_start_pos(bkp_dict['pos1'])
            bkp.set_stop_pos(bkp_dict['pos2'])

        bkp.set_variant_site_1(bkp.get_pos())
        bkp.set_variant_site_2(bkp.get_stop_pos())

        if bkp.get_desc().startswith(EXONIC_IDENTIFIER):
            bkp.set_exon(bkp.get_desc().split(" ")[1])
            bkp.set_site("exon " + str(bkp.get_exon()))
            if __bkp_type != 1:
                bkp.set_variant_site_1(bkp.get_start_pos())
                bkp.set_variant_site_2(bkp.get_pos())

        elif bkp.get_desc().startswith(INTRONIC_IDENTIFIER):
            __exon = int(bkp.get_desc().split(" ")[5])
            bkp.set_exon(__exon)
            bkp.set_intron(__exon)

            if AFTER_IDENTIFIER in bkp.get_desc():
                if __bkp_type == 1:
                    bkp.set_exon(__exon + 1)
                else:
                    bkp.set_variant_site_1(bkp.get_start_pos())
                    bkp.set_variant_site_2(bkp.get_pos())

            if BEFORE_IDENTIFIER in bkp.get_desc():
                bkp.set_intron(__exon - 1)
                if __bkp_type != 1:
                    bkp.set_exon(__exon - 1)
                    bkp.set_variant_site_1(bkp.get_start_pos())
                    bkp.set_variant_site_2(bkp.get_pos())

            bkp.set_site("intron {}".format(bkp.get_intron()))
        elif bkp.get_transcript() == "NM_004449":
            bkp.set_exon(4)
            bkp.set_site("intron 3")
            logger.warning("Breakpoint attributes are estimated for %s:%s. Manual review required" % (bkp.get_chrom(), bkp.get_pos()))
    except Exception as e:
        logger.error("Exception in get_bkp_info")
        logger.error(e)

def get_bkp_type(bkp,fusion,orientation):
    if any([fusion == 1 and orientation == 2,
            fusion == 0 and bkp.get_strand() == POSITIVE_STRAND_IDENTIFIER and orientation == 1,
            fusion == 0 and bkp.get_strand() == NEGATIVE_STRAND_IDENTIFIER and orientation == 2]):
        return 1
    elif any([fusion == 1 and orientation == 1,
            fusion == 0 and bkp.get_strand() == POSITIVE_STRAND_IDENTIFIER and orientation == 2,
            fusion == 0 and bkp.get_strand() == NEGATIVE_STRAND_IDENTIFIER and orientation == 1]):
        return 2
    return 0

def get_bkpsite_note(fusion, __b1, __b2, logger):
    """
    Determine breakpoint details.
    sv, (None or bkp), (None or bkp) -> str
    """
    __note1 = "One of the breakpoints is within"
    __note2 = "The breakpoints are within"
    __note = ""
    try:
        if __b1 is not None and __b2 is not None and __b1.get_site().lower().startswith(EXONIC_IDENTIFIER.lower()) and __b2.get_site().lower().startswith(EXONIC_IDENTIFIER.lower()):
            __note = "{} {} {} and {} {}.".format(__note2, __b1.get_gene(), __b1.get_site(), __b2.get_gene(), __b2.get_site())
            if __b1.get_gene() == __b2.get_gene():
                __note = "{} {}.".format(__note2, __b1.get_site()) if __b1.get_site() == __b2.get_site() else "{} {} and {}.".format(__note2, __b1.get_site(), __b2.get_site())
        else:
            if __b1 is not None and __b1.get_site().lower().startswith(EXONIC_IDENTIFIER.lower()):
                __note = "{} {} {}.".format(__note1, __b1.get_gene(), __b1.get_site()) if fusion else "{} {}.".format(__note1, __b1.get_site())
            elif __b2 is not None and __b2.get_site().lower().startswith(EXONIC_IDENTIFIER.lower()):
                __note = "{} {} {}.".format(__note1, __b2.get_gene(), __b2.get_site()) if fusion else "{} {}.".format(__note1, __b2.get_site())
        return __note
    except Exception as e:
        logger.error("Exception in get_bkpsite_note")
        logger.error(e)
        return  ""

def get_exon_order(bkp, order):
    """
    Determine gene order for notes
    bkp, int -> str
    """
    if any([(bkp.get_strand() == POSITIVE_STRAND_IDENTIFIER and order == 1), (bkp.get_strand() == NEGATIVE_STRAND_IDENTIFIER and order == 2), (order == 4)]): 
        ordert = (bkp.get_exon(), bkp.get_last_exon())
    elif any([(bkp.get_strand() == NEGATIVE_STRAND_IDENTIFIER and order == 1), (bkp.get_strand() == POSITIVE_STRAND_IDENTIFIER and order == 2), (order == 3)]): 
        ordert = (bkp.get_first_exon(), bkp.get_exon())
    else:
        raise Exception("Something terribly went wrong with determining exon order for notes!")
    if ordert[0] == ordert[1]: return "exon %s" % (ordert[0])
    return "exons %s - %s" % ordert

def get_exons_involved(__sv,__bkp_1,__bkp_2,logger):
    """
    Get exons involved in an sv object based on the variant type
    and the breakpoint sites
    sv -> None
    """
    __bkp_1, __bkp_2 = __sv.get_annotation_partner_1(), __sv.get_annotation_partner_2()
    if __sv.get_is_protein_fusion():
        get_bkp_info(__bkp_1, 1, 1, logger)
        get_bkp_info(__bkp_2, 2, 1, logger)
        __sv.set_bkpsites(get_bkpsite_note(__sv.get_is_protein_fusion(),__bkp_1, __bkp_2,logger))
        __note1 = __bkp_1.get_gene() + " " + get_exon_order(__bkp_1, 3)
        __note2 = __bkp_2.get_gene() + " " + get_exon_order(__bkp_2, 4)
        if __sv.get_is_known_oncokb_fusion(): return "%s and %s." % (__note1, __note2)
        return "%s to %s." % (__note1, __note2)
    elif all([__bkp_1.get_is_panel(), __bkp_2.get_is_panel(), __bkp_1.get_is_coding(), __bkp_2.get_is_coding()]):
        get_bkp_info(__bkp_1, 1, 0, logger)
        get_bkp_info(__bkp_2, 2, 0, logger)
        if __sv.check_if_sv_type_translocation():
            return "%s %s and %s %s." % (__bkp_1.get_gene(), __bkp_1.get_site(),__bkp_2.get_gene(), __bkp_2.get_site())
        elif __sv.check_is_intragenic():
            __intra1, __intra2 = (1, 2) if __bkp_1.get_strand() == POSITIVE_STRAND_IDENTIFIER else (2, 1)
            get_bkp_info(__bkp_1, __intra1, 0, logger)
            get_bkp_info(__bkp_2, __intra2, 0, logger)
            __sv.set_bkpsites(get_bkpsite_note(__sv.get_is_protein_fusion(),__bkp_1,__bkp_2,logger))
            if __bkp_1.get_intron() == __bkp_2.get_intron() != None:
                __note1 = "intron %s" % (__bkp_1.get_intron())
            elif not __bkp_1.get_exon() == __bkp_2.get_exon(): 
                __note1 = "exons %s - %s" % (__bkp_1.get_exon(),__bkp_2.get_exon())  
            else: __note1 = "exon %s" % (__bkp_1.get_exon())
            __bkp_1.set_variant_site_1(__bkp_1.get_pos())
            __bkp_2.set_variant_site_1(__bkp_1.get_pos())
            __bkp_1.set_variant_site_2(__bkp_2.get_pos())
            __bkp_2.set_variant_site_2(__bkp_2.get_pos())
            return "%s." % (__note1)
        else:
            __sv.set_bkpsites(get_bkpsite_note(__sv.get_is_protein_fusion(),__bkp_1, __bkp_2,logger))
            __note1 = __bkp_1.get_gene() + " " + get_exon_order(__bkp_1, 1)
            __note2 = __bkp_2.get_gene() + " " + get_exon_order(__bkp_2, 2)
            return "%s and %s." % (__note1, __note2)
    elif __bkp_1.get_is_panel() and __bkp_1.get_is_coding():
        get_bkp_info(__bkp_1, 1, 0, logger)
        if __sv.check_if_sv_type_translocation():
            __note1 = "%s" % (__bkp_1.get_site())
        else:
            __sv.set_bkpsites(get_bkpsite_note(__sv.get_is_protein_fusion(),__bkp_1, None,logger))
            __note1 = get_exon_order(__bkp_1, 1)
        return "%s." % (__note1)
    else:
        get_bkp_info(__bkp_2, 2, 0, logger)
        if __sv.check_if_sv_type_translocation():
            __note1 = "%s" % (__bkp_2.get_site())
        else:
            __sv.set_bkpsites(get_bkpsite_note(__sv.get_is_protein_fusion(), None, __bkp_2,logger))
            __note1 = get_exon_order(__bkp_2, 2)
        return "%s." % (__note1)
    return ""

def get_kinase_status(bkp):
    """
    Get kinase domain annotation for a given bkp
    bkp, df -> (str or None)
    """
    if not bkp.get_is_kinase() : return None
    kinase_interval,interval = [[]]*2
    try:
        interval = [bkp.get_variant_site_1(), bkp.get_variant_site_2()]
        if bkp.get_strand() == NEGATIVE_STRAND_IDENTIFIER:
            interval.sort()
        kinase_interval = KINASE_ANNOTATION[
            (KINASE_ANNOTATION['HUGO'].values == bkp.get_gene())].\
            iloc[0, [1, 4]].values.tolist()
        kinase_interval.sort()
    except AttributeError as e:
        raise Exception(("Breakpoint site positions cannot be determined " \
        "for gene %s. Check if the raw SV call was made on " \
        "Canonical transcript.") % (bkp.gene))
    except (IndexError, NameError, ValueError, Exception) as e:
        pass
    if kinase_interval and interval:
        interval, kinase_interval = \
            map(lambda x: [int(y) for y in x], (interval, kinase_interval))
        overlap = get_overlap(interval, kinase_interval)
        if overlap == 0:
            return "does not include the kinase domain of %s" % (bkp.get_gene())
        else:
            if min(interval) <= min(kinase_interval) and \
                    max(interval) >= max(kinase_interval):
                bkp.set_is_entire_kinase(1)
                return "includes the kinase domain of %s" % (bkp.get_gene())
            else:
                bkp.set_is_entire_kinase(-1)
                return "includes a part of the kinase domain of %s" % (bkp.get_gene())
    return None

def get_misc_notes(sv):
    """
    get frame and kinase domain information for
    fusion variants only
    sv -> str
    """

    __prefix = "The fusion" if sv.get_is_protein_fusion() else "The rearrangement"
    misc_note,frame_note = "",None
    if sv.get_is_protein_fusion() and sv.check_if_sv_in_frame():
        frame_note = "is predicted to be in frame"
    if sv.check_is_intragenic():
        kinase1, kinase2 = get_kinase_status(sv.get_annotation_partner_1()), None
    else:
        kinase1, kinase2 = map(lambda x: get_kinase_status(x), (sv.get_annotation_partner_1(), sv.get_annotation_partner_2()))
    
    if kinase1 and kinase2 and (sv.get_annotation_partner_1().get_is_entire_kinase() == sv.get_annotation_partner_2().get_is_entire_kinase()):  
        regex_str = r'\b' + 'does not include the kinase domain of ' + r'\b|' +\
            r'\b' + 'includes the kinase domain of ' + r'\b|' +\
            r'\b' + 'includes a part of the kinase domain of ' + r'\b'
        kinase2 = re.sub(regex_str, "", kinase2, count=1)
        kinase1 = kinase1.replace("domain", "domains").replace("a part", "parts")

    misc_note = " and ".join(filter(lambda x: isinstance(x, str),[frame_note, kinase1, kinase2]))

    sv.overwrite_protein_fusion()

    if sv.check_if_sv_type_translocation() and not sv.get_is_protein_fusion() : return ""

    if misc_note:
        return "{prefix} {note}.".format(prefix=__prefix, note=misc_note)
    else:
        return ""

def clinical_warning_note(sv):
    bkps_involved = filter(lambda x: x.get_is_kinase_hotspot() ,[sv.get_annotation_partner_1(), sv.get_annotation_partner_2()])
    bkps_genes = list(set(map(lambda x: x.get_gene(), bkps_involved)))
    genes_involved = " and ".join(bkps_genes)
    warning_note = (
        "Functional significance is undetermined and "
        "further testing to determine the presence "
        "or absence of a targetable oncogenic fusion involving %s "
        "is required. This sample has been nominated for further "
        "analysis using the Archer targeted RNAseq assay. "
        "Archer will be performed and reported under a separate "
        "accession number if additional material is available."
    ) % (genes_involved)
    basic_warning_note = (
        "Functional significance is undetermined. "
        "Further testing to characterize the structural variant "
        "involving %s is required."
    ) % (genes_involved)
    return warning_note



def functional_significance(sv):
    """
    Determine if functional significance note is neccessary
    based on variant type
    sv -> str
    """
    if sv.get_is_known_oncokb_fusion():
        return ""
    elif sv.get_annotation_partner_1().get_is_kinase_hotspot() or sv.get_annotation_partner_2().get_is_kinase_hotspot():
        return clinical_warning_note(sv)
    elif sv.check_if_sv_type_deletion() and \
        any([ all([sv.get_annotation_partner_1().get_is_panel_and_coding(),sv.get_annotation_partner_1().get_tumor_suppressor()]), \
            all([sv.get_annotation_partner_2().get_is_panel_and_coding(),sv.get_annotation_partner_2().get_tumor_suppressor()])]):
        return ""  
    else:
        return "Functional significance is undetermined."


def get_sv_oncokb_type(sv):
    r"""
    This function determines the Oncokb special type of the SV

    Parameters
    ----------
    sv : SVAnnotator

    Returns
    -------
    sv_type : string

    """
    if sv.get_annotation().startswith("EGFR (NM_005228) rearrangement:"):
        if re.search(r"exons \d+ - \d+", sv.get_exons_note()):
            se = re.search(r"exons \d+ - \d+", sv.get_exons_note())
            exon1, exon2 = se.group().replace("exons ","").split(" - ")
            if exon1 == "2" and exon2 == "7" and not sv.check_if_sv_type_duplication():
                return SV_SPECICAL_TYPE_VIII
            elif exon1 == "25" and exon2 in ("26", "27", "28" ) and not sv.check_if_sv_type_duplication():
                return SV_SPECICAL_TYPE_CTD
            elif exon1 == "18" and exon2 in ("25", "26") and sv.check_if_sv_type_duplication():
                return SV_SPECICAL_TYPE_KDD
    elif (sv.get_annotation_partner_1().get_transcript() == "NM_004449" and sv.get_annotation_partner_1().get_is_coding()) or (sv.get_annotation_partner_2().get_transcript() == "NM_004449" and sv.get_annotation_partner_2().get_is_coding()):
        return SV_SPECICAL_TYPE_ERG
    elif sv.get_annotation_partner_1().get_gene() == "CDKN2A" or sv.get_annotation_partner_2().get_gene() == "CDKN2A":
        __cdkn2a_tx = filter(lambda x: len(x) > 0, set([sv.get_annotation_partner_1().get_transcript(), sv.get_annotation_partner_2().get_transcript()]))
        if __cdkn2a_tx.__len__() == 2: 
            return SV_SPECICAL_TYPE_CDKN2A_BOTH_ISOFORM
        elif __cdkn2a_tx.__len__() == 1: 
            tx = __cdkn2a_tx.pop()
            if tx == "NM_058195": return SV_SPECICAL_TYPE_CDKN2A_14_ISOFORM
            else: return SV_SPECICAL_TYPE_CDKN2A_16_ISOFORM
        return SV_SPECICAL_TYPE_UNK
    elif sv.get_is_protein_fusion(): return SV_SPECICAL_TYPE_FUSION
    return SV_SPECICAL_TYPE_UNK
    
