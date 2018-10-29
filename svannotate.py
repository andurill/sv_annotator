import os
import sys
import requests
import pandas as pd
import numpy as np
import timeit
from main import constants
from main import models
from main import annotation
from main import notes
from main import position


# Global Variables from constants
target_panel = IMPACTv6_targets
panelKinase = IMPACTv6_kinase_targets
oncoKb = OncoKb_known_fusions


def main(args):
    try:
        svtype, bkp1, bkp2, genes, site1, site2, description, connection = args[0].split(",")
    except ValueError:
        raise IncorrectSVInputArguments()

    variant = sv(svtype, bkp1, bkp2, genes, site1, site2, description, connection)

    try:
        variant.__set_variables(target_panel, panelKinase, oncoKb, refFlat)
    except:
        raise

    try:
        annotation = get_variant_annotation(variant)
    except:
        raise
    
    try:
        note = get_variant_note(variant, annotation)
    except:
        raise

    try:
        position = get_variant_position(note)
    except:
        raise
    
    try:
        note = special_cases(note)
    except:
        raise

    result = {
        'Note' : note,
        'Annotation' : annotation,
        'Position' : position
    }

    return result