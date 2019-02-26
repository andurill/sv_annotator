#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import sys
import requests
import pandas as pd
import numpy as np
import logging
import StringIO
from main.constants import *
from main.models import sv
from main.annotation import get_variant_annotation
from main.notes import get_notes


__author__ = "Gowtham J"
__version__ = "1.0.0"
__maintainer__ = "Gowtham J"
__email__ = "jayakumg@mskcc.org"
__status__ = "Development"


# Global Variables from constants
oncoKb = OncoKb_known_fusions
refFlat = refFlat_canonical
tumourSuppressor = IMPACT_TumourSuppressors
hotspot = IMPACT_Hotspots


def annotate_SV(raw):
    note, annotation, position = [None]*3
    try:
        svtype, bkp1, bkp2, genes, site1, \
            site2, description = raw.split(",")
    except ValueError as e:
        logger.critical(e)
        return note, annotation, position

    try:
        variant = sv(svtype, bkp1, bkp2, genes, site1, site2, description)
    except Exception as e:
        logger.critical(e)
        return note, annotation, position

    try:
        variant.expand(transcript_reference, kinase_annotation,
                         hotspot, tumourSuppressor, oncoKb)
    except Exception as e:
        logger.critical(e)
        return note, annotation, position

    try:
        annotation = get_variant_annotation(variant)
        sv.annotation = annotation
    except Exception as e:
        logger.critical(e)
        return note, annotation, position

    try:
        annotation, note, position = get_notes(variant, refFlat_summary, kinase_annotation)
    except Exception as e:
        logger.critical(e)
        return note, annotation, position

    return note, annotation, position


if __name__ == "__main__":
    # Create the logger
    logger = logging.getLogger('basic_logger')
    logger.setLevel(logging.DEBUG)

    # Setup the console handler with a StringIO object
    log_capture_string = StringIO.StringIO()
    ch = logging.StreamHandler(log_capture_string)
    ch.setLevel(logging.DEBUG)

    # Add the console handler to the logger
    logger.addHandler(ch)

    # Call main function
    note, annotation, position = annotate_SV(sys.argv[1])

    # Send log contents to a string and close the stream
    log_contents = log_capture_string.getvalue()
    log_capture_string.close()

    if note and position and annotation:
        status = "SUCCESS"
    else:
        status = "FAILED"
    result = {
        'Note': note,
        'Annotation': annotation,
        'Position': position,
        'status': status,
        'message': log_contents.replace("\n","; ")
    }
    print result
