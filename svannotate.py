#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@Description: Generate human-readable interpretations of structural variants.
@author: jayakumg@mkscc.org
"""


import os
import sys
import requests
import pandas as pd
import numpy as np
import logging
import StringIO
import argparse
import multiprocessing as mp
import contextlib
from main.constants import *
from main.models import sv, build_cache
from main.annotation import get_variant_annotation
from main.notes import get_notes
from main.models import timestamp


# Global Variables from constants
oncoKb = OncoKb_known_fusions
refFlat = refFlat_canonical
tumourSuppressor = IMPACT_TumourSuppressors
hotspot = IMPACT_Hotspots
cache = None


def main():
    parser = argparse.ArgumentParser(description="SV annotator")
    analysis_type = parser.add_mutually_exclusive_group(required=True)
    # single sv processing
    analysis_type.add_argument(
        "-sv",
        "--structural_variant",
        type=str,
        action="store",
        help="Structural variant string",
    )
    # batch processing
    analysis_type.add_argument(
        "-i",
        "--input_file",
        type=argparse.FileType("r"),
        action="store",
        help="input file of multiple structural variants",
    )
    parser.add_argument(
        "-o",
        "--out_file",
        type=argparse.FileType("w"),
        action="store",
        default=os.path.join(os.getcwd(), "SVannotated_results.txt"),
        help="output file to print annotation results for multiple structural variants",
    )
    parser.add_argument(
        "--offline", action="store_true", help="run annotation using cache"
    )
    args = parser.parse_args()

    # Create the logger
    logger = logging.getLogger("basic_logger")
    logger.setLevel(logging.DEBUG)

    # Setup the console handler with a StringIO object
    log_capture_string = StringIO.StringIO()
    ch = logging.StreamHandler(log_capture_string)
    ch.setLevel(logging.DEBUG)

    # Add the console handler to the logger
    logger.addHandler(ch)

    # Call main function
    if args.input_file:
        svdata = pd.read_csv(
            args.input_file,
            header="infer",
            sep="\t",
            dtype=str,
            usecols=[
                "TumorId",
                "NormalId",
                "Chr1",
                "Pos1",
                "Chr2",
                "Pos2",
                "SV_Type",
                "Gene1",
                "Gene2",
                "Site1Description",
                "Site2Description",
                "Fusion",
            ],
        )
        cache_input_data = pd.concat(
            [
                svdata[["Chr1", "Pos1"]].rename(
                    columns={"Chr1": "#CHROM", "Pos1": "POS"}
                ),
                svdata[["Chr2", "Pos2"]].rename(
                    columns={"Chr2": "#CHROM", "Pos2": "POS"}
                ),
            ],
            ignore_index=True,
            axis=0,
        )
        cache_input_data.drop_duplicates(inplace=True)
        cache_input_data["ID"], cache_input_data["REF"], cache_input_data[
            "ALT"
        ] = ".,N,-".split(",")

        try:
            print(timestamp() + "Proceeding to build annotation cache...")
            global cache
            cache = build_cache(
                cache_input_data[["#CHROM", "POS", "ID", "REF", "ALT"]],
                transcript_reference,
            )
        except Exception as e:
            print(timestamp() + "Failed to build annotation cache using VEP.")
            raise

        svdata["coord1"] = svdata["Chr1"] + ":" + svdata["Pos1"]
        svdata["coord2"] = svdata["Chr2"] + ":" + svdata["Pos2"]
        svdata["Genes"] = svdata["Gene1"] + " / " + svdata["Gene2"]
        svdata["SV_Type"] = svdata["SV_Type"].apply(
            lambda x: "INVERSION"
            if x == "INV"
            else (
                "TRANSLOCATION"
                if x == "TRA"
                else ("DELETION" if x == "DEL" else "DUPLICATION")
            )
        )
        svdata = svdata.drop(["Chr1", "Chr2", "Pos1", "Pos2"], axis=1)
        col = [
            u"SV_Type",
            u"coord1",
            u"coord2",
            u"Genes",
            u"Site1Description",
            u"Site2Description",
            u"Fusion",
        ]
        sv_strings = (
            svdata[col]
            .apply(lambda row: ",".join(row.values.astype(str)), axis=1)
            .tolist()
        )

        cores = int(mp.cpu_count() - 1)
        try:
            # create jobs
            print(timestamp() + "Starting variant annotation...")
            with contextlib.closing(mp.Pool(processes=cores)) as pool:
                annotated_SVs = pool.map(annotate_SV, sv_strings)
        except Exception as e:
            print(e)

        new = pd.concat(
            [
                svdata,
                pd.DataFrame(annotated_SVs, columns=["Note", "Annotation", "Position"]),
            ],
            axis=1,
        )

        new.to_csv(args.out_file, header=True, sep="\t", index=False)
        print(timestamp() + "Completed variant annotation!")
    # note, annotation, position = annotate_SV(args.sv, logger)

    # Send log contents to a string and close the stream
    # log_contents = log_capture_string.getvalue()
    # log_capture_string.close()

    # if note and position and annotation:
    #     status = "SUCCESS"
    # else:
    #     status = "FAILED"
    # result = {
    #     "Note": note,
    #     "Annotation": annotation,
    #     "Position": position,
    #     "status": status,
    #     "message": log_contents.replace("\n", "; "),
    # }
    # print(result)


def annotate_SV(raw):
    """
    Main function to initialize sv and breakpoints
    objects based on given inputs and call methods to 
    generate annotation and notes.
    str -> tuple
    """
    note, annotation, position, sv_type = [None] * 4
    try:
        svtype, bkp1, bkp2, genes, site1, site2, description = raw.split(",")
    except ValueError as e:
        return note, annotation, position, sv_type

    try:
        variant = sv(svtype, bkp1, bkp2, genes, site1, site2, description
        variant.expand(
            transcript_reference,
            kinase_annotation,
            hotspot,
            tumourSuppressor,
            oncoKb,
            cache,
        )
        
        annotation = get_variant_annotation(variant)
        
        note, position, sv_type = get_notes(
            variant, refFlat_summary, kinase_annotation
        )
    except Exception as e:
        return note, annotation, position, sv_type

    return (note, annotation, position, sv_type)


if __name__ == "__main__":
    main()
