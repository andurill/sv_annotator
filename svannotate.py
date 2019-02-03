#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os, sys, requests
import pandas as pd
import numpy as np
import timeit
from main.constants import *
from main.models import *
#import main.notes as notes
#import main.position as pos
from main.annotation import *


__author__ = "Gowtham J"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Gowtham J"
__email__ = "jayakumg@mskcc.org"
__status__ = "Development"


# Global Variables from constants
target_panel = IMPACTv6_targets
panelKinase = IMPACTv6_kinase_targets
oncoKb = OncoKb_known_fusions
refFlat = refFlat_canonical
cb_df = cb_df
tumorsuppressor = IMPACT_TumourSuppressors
hotspot = IMPACT_KinaseHotspots


def annotate_SV(raw):
    message, note, annotation, position = [None]*4
    try:
        svtype, bkp1, bkp2, genes, site1, site2, description, connection = raw.split(",")
    except ValueError as e:
        message = e
        return message, note, annotation, position

    try:
        variant = sv(svtype, bkp1, bkp2, genes, site1, site2, description, connection)
    except Exception as e:
        message = e
        return message, note, annotation, position       

    try:
        variant.expand(target_panel, panelKinase, oncoKb, refFlat, hotspot, tumorsuppressor)
    except Exception as e:
        message = e
        return message, note, annotation, position

    try:
        annotation = get_variant_annotation(variant)
        sv.annotation = annotation
    except Exception as e:
        message = e
        return message, note, annotation, position

    return message, note, annotation, position


if __name__ == "__main__":
    message, note, annotation, position = annotate_SV(sys.argv[1])
    if note and position and annotation:
        status = "SUCCESS"
    else:
        status = "FAILED"
    result = {
        'Note': note,
        'Annotation': annotation,
        'Position': position,
        'status': status,
        'message': message
    }
    print result
