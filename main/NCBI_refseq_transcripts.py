#!/usr/bin/env python

"""
@Description: Generate list of all refseq and canonical refseq transcripts from NCBI refseq transcript data in gff3 format
"""

import os, sys, re
from collections import OrderedDict
from operator import itemgetter

patterns = OrderedDict(
    [
        ("GeneID", re.compile(r".*GeneID:([0-9]+)")),
        ("HGNC", re.compile(r".*HGNC:([0-9]+)")),
        ("Gene", re.compile(r".*gene=(\S+);")),
        ("Tx", re.compile(r".*transcript_id=(NM_[0-9]+\.[0-9]+)")),
        ("Tag", re.compile(r".*tag=(.+);")),
    ]
)


def get_IDs(val):
    matches = ["", "", "", "", ""]
    for i, key in enumerate(patterns.keys()):
        try:
            matches[i] = patterns[key].match(val).groups()[0]
        except AttributeError:
            matches[i] = ""
    return matches


refseq_source = sys.argv[1]
with open(refseq_source, "r") as f:
    data = f.read().splitlines()


# select for mRNA entries
data = list(
    filter(
        lambda x: not (x.startswith("#"))
        and x.split("\t")[1] == "BestRefSeq"
        and x.split("\t")[2] == "mRNA",
        data,
    )
)
# reformat
data = [
    "\t".join(list(itemgetter(*[0, 2, 6])(x.split("\t"))) + get_IDs(x.split("\t")[8]))
    for x in data
]
# add header
header = """#!processor NCBI annotwriter
#!genome-build GRCh37.p13
#!genome-build-accession NCBI_Assembly:GCF_000001405.25
#!annotation-date 09/06/2019
#!annotation-source NCBI Homo sapiens Updated Annotation Release 105.20190906
##sequence-region NC_000001.10 1 249250621
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
##source ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
##source https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml
##reference for canonical (Refseq Select) https://www.ncbi.nlm.nih.gov/refseq/refseq_select/
##reference for canonical (MANE Select) https://www.ncbi.nlm.nih.gov/refseq/MANE/
##Sequence_Region   Identifier  Strand  Gene_ID  HGNC    Gene_Name   Transcript_ID  Isoform_Type"""

with open("GRCh37_genomic_release_v105.20190906.txt", "w") as f:
    f.write("\n".join([header] + data))

# select canonical
data = list(filter(lambda x: x.endswith("Select"), data))
with open("GRCh37_genomic_release_v105.20190906_canonical.txt", "w") as f:
    f.write("\n".join([header] + data))
