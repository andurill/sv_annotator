import os, sys
from constants import cb_df


def get_variant_annotation(sv):
    if sv.svtype == "TRANSLOCATION":
        return get_translocation(sv)
    else:
        return get_other_svs(sv)


def get_translocation(sv):
    fusion_type = "rearrangement"
    Annotation = ""
    if sv.isFusion and sv.bkp1.isCoding and sv.bkp2.isCoding:
        if sv.isKnownFusion is True:
            fusion_type = "fusion"
        Annotation = "%s (%s) - %s (%s) %s: " %\
            (sv.fusionPartner1.gene, sv.fusionPartner1.transcript,
             sv.fusionPartner2.gene, sv.fusionPartner2.transcript, fusion_type)
    elif sv.bkp1.isPanel and sv.bkp2.isPanel and sv.bkp1.isCoding and sv.bkp2.isCoding:
        Annotation = "%s (%s) - %s (%s) %s: " %\
            (sv.annotationPartner1.gene, sv.annotationPartner1.transcript,
             sv.annotationPartner2.gene, sv.annotationPartner2.transcript, fusion_type)
    elif sv.bkp1.isPanel and sv.bkp1.isCoding:
        Annotation = "%s (%s) %s: " %\
            (sv.annotationPartner1.gene, sv.annotationPartner1.transcript, fusion_type)
    else:
        Annotation = "%s (%s) %s: " %\
            (sv.annotationPartner2.gene, sv.annotationPartner2.transcript, fusion_type)

    cband1 = get_cytoband(sv.annotationPartner1)
    cband2 = get_cytoband(sv.annotationPartner2)
    t_format = sv.annotationPartner1.chrom, sv.annotationPartner2.chrom,\
        cband1, cband2, sv.annotationPartner1.chrom, sv.annotationPartner1.pos,\
        sv.annotationPartner2.chrom, sv.annotationPartner2.pos
    coordinate = "t(%s,%s)(%s;%s)(chr%s:g.%s::chr%s:g.%s)" % t_format
    return Annotation + coordinate


def get_other_svs(sv):
    svtype = reformat(sv.svtype)
    fusion_type = "rearrangement"
    if sv.isFusion and sv.bkp1.isCoding and sv.bkp2.isCoding:
        if sv.isKnownFusion:
            fusion_type = "fusion"
        gene1, tx1, cdna1 = sv.annotationPartner1.gene, \
            sv.annotationPartner1.transcript, sv.annotationPartner1.cdna
        gene2, tx2, cdna2 = sv.annotationPartner2.gene, \
            sv.annotationPartner2.transcript, sv.annotationPartner2.cdna
        Annotation = "%s (%s) - %s (%s) %s: %s:%s_%s:%s%s" %\
            (gene1, tx1, gene2, tx2, fusion_type,
             cdna1, gene1, cdna2, gene2, svtype)
        print Annotation
        return Annotation
    elif sv.bkp1.isPanel and sv.bkp2.isPanel and sv.bkp1.isCoding and sv.bkp2.isCoding:
        gene1, tx1, cdna1 = sv.bkp1.gene, sv.bkp1.transcript, sv.bkp1.cdna
        gene2, tx2, cdna2 = sv.bkp2.gene, sv.bkp2.transcript, sv.bkp2.cdna
        if sv.isIntragenic:
            Annotation = "%s (%s) %s: %s_%s:%s" %\
                (gene1, tx1, fusion_type, cdna1, cdna2, svtype)
        else:
            Annotation = "%s (%s) - %s (%s) %s: %s:%s_%s:%s%s" %\
                (gene1, tx1, gene2, tx2, fusion_type,
                 cdna1, gene1, cdna2, gene2, svtype)
            print Annotation
            return Annotation
    elif sv.bkp1.isPanel and sv.bkp1.isCoding:
        gene1, tx1, cdna1 = sv.bkp1.gene, sv.bkp1.transcript, sv.bkp1.cdna
        if not cdna1.startswith("chr"):
            cdna2 = "chr" + sv.bkp2.chrom + ":g." + str(sv.bkp2.pos)
            Annotation = "%s (%s) %s: %s:%s_%s%s" %\
                (gene1, tx1, fusion_type, cdna1, gene1, cdna2, svtype)
            print Annotation
            return Annotation
        else:
            raise BreakPointIntergenic(sv.bkp1)
    else:
        gene2, tx2, cdna2 = sv.bkp2.gene, sv.bkp2.transcript, sv.bkp2.cdna
        if not cdna2.startswith("chr"):
            cdna1 = "chr" + sv.bkp1.chrom + ":g." + str(sv.bkp1.pos)
            Annotation = "%s (%s) %s: %s:%s_%s%s" %\
                (gene2, tx2, fusion_type, cdna2, gene2, cdna1, svtype)
            print Annotation
            return Annotation
        else:
            raise BreakPointIntergenic(sv.bkp2)
    return


def get_cytoband(bkp):
    chrom, coord = bkp.chrom, bkp.pos
    which_cytoband = cb_df[(cb_df['Chr'].values == chrom) &
                           (cb_df['Bp_start'].values <= coord) &
                           (cb_df['Bp_stop'].values >= coord)]
    if len(which_cytoband) == 0:
        raise MissingCytoBand(bkp)
    elif len(which_cytoband) > 1:
        raise MultipleCytoBand(bkp)
    locus = which_cytoband.values[0]
    return locus[0] + locus[1]


def reformat(svtype):
    svdict = {
        "DUPLICATION": "dup",
        "DELETION": "del",
        "INVERSION": "inv"
    }
    return svdict[svtype]
