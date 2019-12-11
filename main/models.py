#!/usr/bin/env python2
import os
import sys
import logging
import requests
import pandas as pd
from config import VEP, VEP_CACHE, PERL, FASTA
import random
import subprocess
from datetime import datetime

logger = logging.getLogger("basic_logger")


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S ")


def create_cache(input_data):
    input_data["ID"], input_data["REF"], input_data["ALT"] = zip(
        *[(".", "N", "-")] * len(input_data.index)
    )
    input_data[["#CHROM", "POS", "ID", "REF", "ALT"]].to_csv(
        "tmp.vcf", header=True, sep="\t", index=False
    )

    return


class bkp(object):
    """
    Class to represent a breakpoint and other related attributes and features.
    """

    def __init__(self, chrom, pos, gene, desc):
        try:
            self.chrom = str(chrom)
            self.pos = int(pos)
            self.gene = str(gene)
            self.desc = str(desc)
        except TypeError:
            raise Exception(
                "Could not create a new instance of bkp class due to inappropriate values for parameters."
            )
        except Exception:
            raise Exception("Unexpected error:", sys.exc_info()[0])

    def expand(
        self, transcript_reference, kinase_annotation, hotspot, tumourSuppressor, cache
    ):
        reference = transcript_reference[
            transcript_reference["Gene"] == self.gene
        ].to_dict("list")
        self.transcript = reference["Lookup_Transcript"] or None
        self.cdna = "chr" + self.chrom + ":g." + str(self.pos)
        if self.transcript is None:
            logger.warning("Cannot find canonical transcript for " + str(self.gene))
            # self.cdna = "chr" + self.chrom + ":g." + str(self.pos)
        else:
            try:
                self.transcript, self.cdna = get_cdna_pos(self, cache)
                # self.transcript = self.transcript.split(".")[0]
                # if self.transcript == "ENST00000331340":
                #    self.transcript = "NM_006060"
                self.transcript = dict(
                    zip(
                        reference["Lookup_Transcript"], reference["Reported_Transcript"]
                    )
                )[self.transcript]
            except Exception as e:
                logger.warning(e)

        if "(+)" in self.desc:
            self.strand = "+"
        elif "(-)" in self.desc:
            self.strand = "-"
        else:
            raise Exception(
                "Cannot get strand orientation info for breakpoint %s:%s."
                % (self.chrom, self.pos)
            )

        if self.cdna and self.cdna.startswith("c."):
            self.isCoding = True
        else:
            self.isCoding = False

        self.isPanel = False
        try:
            if set(reference["IsPanel"]).pop() == 1:
                self.isPanel = True
        except KeyError:
            logger.warning(
                "%s not present in the reference canonical transcript file."
                % (self.gene)
            )

        if all([self.gene in tumourSuppressor, self.isPanel, self.isCoding]):
            self.isTumourSuppressor = True
        else:
            self.isTumourSuppressor = False

        if self.gene in hotspot and self.isCoding:
            self.isHotspot = True
        else:
            self.isHotspot = False

        if self.gene in kinase_annotation["HUGO"].tolist() and self.isCoding:
            self.isKinase = True
        else:
            self.isKinase = False


class sv(object):
    """
    Class to represent a structural variant and other related attributes and features.
    """

    message = ""

    def __init__(self, svtype, bkp1, bkp2, genes, site1, site2, description):
        self.svtype = svtype
        self.site1 = site1
        self.site2 = site2
        self.description = description
        try:
            self.chr1, self.pos1 = bkp1.split(":")  # IncorrectBkpFormat
            self.chr2, self.pos2 = bkp2.split(":")  # IncorrectBkpFormat
            self.gene1, self.gene2 = genes.split(" / ")  # IncorrectGenesFormat
        except ValueError:
            raise Exception(
                "Could not create a new instance of sv class due to incorrect format of arguments."
            )

    def expand(
        self,
        transcript_reference,
        kinase_annotation,
        hotspot,
        tumourSuppressor,
        oncoKb,
        cache,
    ):
        self.bkp1 = bkp(self.chr1, self.pos1, self.gene1, self.site1)
        self.bkp2 = bkp(self.chr2, self.pos2, self.gene2, self.site2)
        self.bkp1.expand(
            transcript_reference, kinase_annotation, hotspot, tumourSuppressor, cache
        )
        self.bkp2.expand(
            transcript_reference, kinase_annotation, hotspot, tumourSuppressor, cache
        )

        # Define key variables
        if not (self.bkp1.isPanel or self.bkp2.isPanel):
            raise GenesNotInPanel()

        if (self.bkp1.isPanel and self.bkp1.isCoding) or (
            self.bkp2.isPanel and self.bkp2.isCoding
        ):
            pass
        else:
            raise BothBreakpointsNoncoding()
        # Is SV intragenic?
        if all(
            [
                self.bkp1.gene == self.bkp2.gene,
                self.bkp1.transcript == self.bkp2.transcript,
                self.bkp1.isCoding,
                self.bkp2.isCoding,
            ]
        ):
            self.isIntragenic = True
        else:
            self.isIntragenic = False

        # Fusion variables
        if (
            (
                self.description.startswith("Protein Fusion")
                or self.description.startswith("Transcript Fusion")
            )
            and self.bkp1.isCoding
            and self.bkp2.isCoding
        ):
            self.isFusion = True
            s = self.description
            self.fusionGene = s[s.find("{") + 1 : s.find("}")]
            # check if fusion is defined in the correct format Gene1:Gene2
            try:
                self.fusionPartner1, self.fusionPartner2 = self.fusionGene.split(":")
            except ValueError:
                raise IncorrectDescriptionFormat()

            # Check if fusion partners match the genes provided as inputs
            if (
                self.fusionPartner1 == self.bkp1.gene
                and self.fusionPartner2 == self.bkp2.gene
            ):
                self.fusionPartner1, self.fusionPartner2 = self.bkp1, self.bkp2
            elif (
                self.fusionPartner2 == self.bkp1.gene
                and self.fusionPartner1 == self.bkp2.gene
            ):
                self.fusionPartner2, self.fusionPartner1 = self.bkp1, self.bkp2
            else:
                raise FusionGeneConflict()

            # Is fusion known?
            if self.fusionGene in oncoKb:
                self.isKnownFusion = True
            else:
                self.isKnownFusion = False
        else:
            self.fusionGene = None
            self.isFusion = False
            self.isKnownFusion = False
            if "Fusion:" in self.description:
                if not self.bkp1.isCoding:
                    message = "coding cDNA annotation cannot be found for %s %s." % (
                        self.bkp1.gene,
                        self.bkp1.transcript,
                    )
                elif not self.bkp2.isCoding:
                    message = "coding cDNA annotation cannot be found for %s %s." % (
                        self.bkp2.gene,
                        self.bkp2.transcript,
                    )
                logger.warning(message)

        # Annotation variables
        self.annotationPartner1, self.annotationPartner2 = self.bkp1, self.bkp2
        if any(
            [
                self.svtype == "TRANSLOCATION"
                and any(
                    [
                        self.bkp1.chrom.isdigit()
                        and self.bkp2.chrom.isdigit()
                        and int(self.bkp1.chrom) > int(self.bkp2.chrom),
                        str(self.bkp1.chrom) == "X",
                        str(self.bkp1.chrom) == "Y",
                    ]
                ),
                self.isIntragenic and self.bkp1.strand == "-",
            ]
        ):
            self.annotationPartner1, self.annotationPartner2 = self.bkp2, self.bkp1
        elif all([self.isFusion, self.bkp1.isCoding, self.bkp2.isCoding]):
            self.annotationPartner1, self.annotationPartner2 = (
                self.fusionPartner1,
                self.fusionPartner2,
            )


def get_cdna_pos(bkp, cache):
    """
    Get cdna position for a bkp object by querying
    vep server
    bkp -> tuple
    """
    if cache is None:
        tx, cdna = get_cdna_pos_from_api(bkp)
    else:
        tx, cdna = get_cdna_pos_from_cache(bkp, cache)
    return tx, cdna
    # return (
    #     get_cdna_pos_from_cache(bkp, cache)
    #     if cache is not None
    #     else get_cdna_pos_from_api(bkp)
    # )


def build_cache(bkps, transcript_reference):
    # print(transcript_reference["Lookup_Transcript"])

    def filter_for_select_tx(list_of_dict):
        try:
            return list(
                filter(
                    lambda x: x["transcript_id"]
                    in transcript_reference["Lookup_Transcript"].values,
                    list_of_dict,
                )
            )
        except TypeError as e:
            return []

    try:
        tmp_dir = os.environ["TMP"]
    except KeyError:
        tmp_dir = os.getcwd()
    input_file_name = "_".join(["variant", str(random.getrandbits(128))]) + ".vcf"
    out_file_name = input_file_name.replace(".vcf", ".txt")
    with open(os.path.join(tmp_dir, input_file_name), "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        bkps.to_csv(f, header=True, index=None, sep="\t", mode="a")

    vep_cmd = " ".join(
        [
            PERL,
            VEP,
            "--species",
            "homo_sapiens_merged",
            "--assembly",
            "GRCh37",
            "--offline",
            "--no_progress",
            "--no_stats",
            "--buffer_size",
            "10000",
            "--hgvs",
            "--minimal",
            "--canonical",
            "--dir",
            VEP_CACHE,
            "--fasta",
            FASTA,
            "--fork",
            "12",
            "--json",
            "--input_file",
            os.path.join(tmp_dir, input_file_name),
            "--out",
            os.path.join(os.getcwd(), out_file_name),
        ]
    )

    try:
        print(timestamp() + "Running VEP...")
        subprocess.check_call(vep_cmd, shell=True)
    except subprocess.CalledProcessError:
        raise
    try:
        print(timestamp() + "Reading json outputs from vep:")
        annotation_results = pd.read_json(
            os.path.join(os.getcwd(), out_file_name), lines=True
        )[["id", "transcript_consequences"]]
    except IOError as e:
        print(timestamp() + "VEP annotated json file cannot be found. Check it VEP ran successfully.")
        raise

    print(timestamp() + "Parsing JSON results.")
    annotation_results.dropna(inplace=True)
    annotation_results["id"] = (
        annotation_results["id"].str.replace("_N/-", "").str.replace("_", ":")
    )
    print(timestamp() + "Filtering for canonical transcripts.")
    annotation_results["transcript_consequences"] = annotation_results[
        "transcript_consequences"
    ].apply(filter_for_select_tx)
    # print(annotation_results["transcript_consequences"])
    print(timestamp() + "Completed filtering")
    annotation_results = (
        pd.DataFrame(
            annotation_results["transcript_consequences"].tolist(),
            index=annotation_results["id"],
        )
        .stack()
        .reset_index(name="transcript_consequences")[["id", "transcript_consequences"]]
    )
    print(timestamp() + "Coverted to long.")
    # print(annotation_results["transcript_consequences"])
    # annotation_results[["hgvsc", "transcript_id", "gene_symbol"]] = annotation_results[
    #     "transcript_consequences"
    # ].apply(pd.Series)[["hgvsc", "transcript_id", "gene_symbol"]]
    annotation_results[["hgvsc", "transcript_id"]] = annotation_results[
        "transcript_consequences"
    ].apply(pd.Series)[["hgvsc", "transcript_id"]]
    print(timestamp() + "Expanded columns.")
    annotation_results.drop(["transcript_consequences"], axis=1, inplace=True)
    # annotation_results.dropna(subset=["gene_symbol"], inplace=True)
    annotation_results["hgvsc"] = (
        annotation_results["hgvsc"].str.replace(r".*:", "").str.replace(r"del.*", "")
    )
    print(timestamp() + "Completed parsing JSON results")
    return annotation_results


def get_cdna_pos_from_cache(bkp, cache):
    cdna = None
    for tx in bkp.transcript:
        try:
            query_cdna = cache[
                (cache["id"].values == ":".join([bkp.chrom, str(bkp.pos)]))
                & (cache["transcript_id"].values == tx)
            ]["hgvsc"].tolist()[0]
            if query_cdna.startswith("c.-") or query_cdna.startswith("c.*"):
                continue
            else:
                cdna = str(query_cdna)
                break
        except (KeyError, IndexError) as e:
            cdna = "chr" + bkp.chrom + ":g." + str(bkp.pos)
    if cdna is None:
        cdna = "chr" + bkp.chrom + ":g." + str(bkp.pos)
    return tx, cdna


def get_cdna_pos_from_api(bkp):
    dummy_ref = "C"
    query = make_query(bkp, dummy_ref)
    cdna = None

    # get request max twice to VEP for annotation
    request = make_get_request(query)
    if not request.ok:
        s = request.text
        actual_ref = s[s.find("(") + 1 : s.find(")")]
        if actual_ref not in ["A", "C", "G", "T"] or actual_ref == dummy_ref:
            request.raise_for_status()
        else:
            query = make_query(bkp, actual_ref)
            request = make_get_request(query)
            if not request.ok:
                request.raise_for_status()
    decoded = request.json()
    try:
        result = dict(str(s).split(":", 1) for s in decoded[0]["hgvsc"])
    except KeyError:
        result = {}
        logger.warning("Cannot find any cDNA annotations for gene " + bkp.gene)
    for tx in bkp.transcript:
        if tx in result:
            query_cdna = result[tx]
            if query_cdna.startswith("c.-") or query_cdna.startswith("c.*"):
                pass
            else:
                cdna = query_cdna[:-3]
                break
    if not cdna:
        cdna = "chr" + bkp.chrom + ":g." + str(bkp.pos)
    return tx, cdna


def make_get_request(query):
    """
    Generate query string in the format required by vep
    str -> request
    """
    server = "http://grch37.rest.ensembl.org"
    ext = "/variant_recoder/human/" + query
    try:
        request = requests.get(
            server + ext, headers={"Content-Type": "application/json"}, timeout=5
        )
    except requests.ConnectionError as e:
        raise e
    except requests.exceptions.RequestException as e:
        raise e
    return request


def make_query(bkp, dummy_ref):
    """
    Make dummy query string for the second
    request sent to vep server
    bkp, str -> str
    """
    chrom, coord = bkp.chrom, bkp.pos
    revcomp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    dummy_alt = revcomp[dummy_ref]
    query = chrom + ":g." + str(coord) + dummy_ref + ">" + dummy_alt + "?"
    return query


class Error(Exception):
    """Base class for other exceptions"""

    pass


class MissingCytoBand(Error):
    """Raised when no cytoband was identified for a breakpoint"""

    def __init__(self, bkp):
        Exception.__init__(self, "No cytobands identified for the breakpoint: " + bkp)


class MultipleCytoBand(Error):
    """Raised when multiple cytoband were identified for a breakpoint"""

    def __init__(self):
        Exception.__init__(
            self, "Multiple cytobands identified for the breakpoint: " + bkp
        )


class CanonicalTranscriptNotFound(Error):
    """Raised when canonical transcript not found"""

    def __init__(self, gene):
        Exception.__init__(
            self, "Cannot find a canonical transcript for the gene: " + gene
        )


class IncorrectGenesFormat(Error):
    """Raised when gene1 / gene2 format does not match expected"""

    def __init__(self):
        Exception.__init__(
            self,
            "Input genes are not in required format. Should be in 'Gene1 / Gene2' format",
        )


class IncorrectBkpFormat(Error):
    """Raised when a bkp format does not match expected"""

    def __init__(self, bkp):
        Exception.__init__(
            self,
            "Format of breakpoint "
            + str(bkp)
            + "does not match the required format. Example '7:1234567'",
        )


class IncorrectDescriptionFormat(Error):
    """Raised when the input fusion description does not match expected"""

    def __init__(self):
        Exception.__init__(
            self,
            "Format of fusion description does not meet the required format. Example 'Protein Fusion: {FGFR3:TACC3}'",
        )


class GenesNotInPanel(Error):
    """Raised when both genes not in panel"""

    def __init__(self):
        Exception.__init__(self, "Neither of the gene partners is in the target panel")


class FusionGeneConflict(Error):
    """Raised when one of the fusion parter does not match either of the input genes"""

    def __init__(self):
        Exception.__init__(
            self, "The fusion partner genes do not match the two breakpoint genes"
        )


class cdnaNotFound(Error):
    """Raised when vep cannot generate cdna annotation for a breakpoint"""

    def __init__(self, bkp):
        Exception.__init__(
            self,
            "Cannot determine cDNA annotation using vep for the breakpoint: " + bkp,
        )


class IncorrectSVInputArguments(Error):
    """Raised when the input arguments to sv class in incorrect"""

    def __init__(self):
        Exception.__init__(
            self,
            "The number arguments provided for sv class in correct. It should be a \
            comma-separater string of 'svtype, bkp1, bkp2, genes, site1, site2, description, connection'",
        )


class BreakPointIntergenic(Error):
    """Raised when breakpoint is in an intergenic region"""

    def __init__(self, bkp):
        Exception.__init__(
            self,
            "Breakpoint is is the intergenic region for the breakpont %s:%s"
            % (bkp.chrom, str(bkp.pos)),
        )


class GenomicPosition(Error):
    """Raised when breakpoint position cannot be determined in the form of genomic or DNA coordinates"""

    def __init__(self, bkp):
        Exception.__init__(
            self,
            "Genomic position cannot be determined in the form of genomic or DNA for breakpoint %s:%s"
            % (bkp.chrom, str(bkp.pos)),
        )


class BothBreakpointsNoncoding(Error):
    """Raised when not even one of the breakpoints is in panel and is in coding region"""

    def __init__(self):
        Exception.__init__(
            self,
            "Atleast one of the breakpoints need to be in target panel and in coding region",
        )

