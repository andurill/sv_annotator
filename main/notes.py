import os, sys
import pandas as pd
import numpy as np
from notes_constants import *


refFlat_summary = pd.read_csv("./data/refFlat_summary.txt", sep="\t", dtype={
                              'a': str, 'b': str, 'c': str, 'd': str, 
                              'e': int, 'f': int, 'g': int})



def get_note(sv):
    note = None
    return note


def get_prefix(sv):
    prefix = "Note: The "
    fusion_type = "rearrangement"
    if sv.isFusion:
        if sv.isKnownFusion:
            fusion_type = "fusion"
        prefix += "%s - %s %"
    return


def get_kinase(sv):
    for bkp in (sv.bkp1, sv.bkp2):
        if bkp.isKinase:

    
