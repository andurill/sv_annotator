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
        raise SVInputArgumentsError()

    variant = SV(svtype, bkp1, bkp2, genes, site1, site2, description, connection)

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


class sv(object):
    def __init__(self, svtype, bkp1, bkp2, genes, site1, site2, description, connection):
       self.__svtype = svtype
       self.__bkp1 = bkp1
       self.__bkp2 = bkp2
       self.__genes = genes
       self.__site1 = site1
       self.__site2 = site2
       self.__description = description
       self.__connection = connection

    
    # Define key variables
    def __set_variables(self, target_panel, panelKinase, oncoKb, refFlat):
        try:
            self.__gene1, self.__gene2 = genes.split(" / ")
        except ValueError:
            raise IncorrectGenesFormat()
        
        try:
            self.__chr1, self.__pos1 = self.__bkp1.split(":")
        except ValueError:
            raise IncorrectBkpFormat()

        try:
            self.__chr2, self.__pos2 = self.__bkp2.split(":")
        except ValueError:
            raise IncorrectBkpFormat()

        # get canonical transcript(s) for gene1, otherwise raise exception
        if refFlat[self.__gene1]:
            self.__transcript1 = refFlat[self.__gene1]
        else:
            raise CanonicalTranscriptNotFound(self.__gene1)

        # get canonical transcript(s) for gene2, otherwise raise exception
        if refFlat[self.__gene2]:
            self.__transcript2 = refFlat[self.__gene2]
        else:
            raise CanonicalTranscriptNotFound(self.__gene2)

        # Is gene1 in panel?
        if target_panel[self.__gene1]:
            self.__isGene1InPanel = True
        else:
            self.__isGene1InPanel = False

        # Is gene2 in panel?
        if target_panel[self.__gene2]:
            self.__isGene2InPanel = True
        else:
            self.__isGene2InPanel = False

        if self.__isGene1InPanel is False and self.__isGene2InPanel is False:
            raise GenesNotInPanel()

        # Is SV intragenic?
        if self.__gene1 == self.__gene2:
            self.__isIntragenic = True
        else:
            self.__isIntragenic = False

        # Are one or both genes in panel and kinase?
        self.__kinase_genes = set()
        for gene in self.__gene1, self.__gene2:
            if target_panel[gene] and panelKinase[gene]:
                self.__kinase_genes.add(gene)

        # IF PROTEIN FUSION
        if self.__description.startswith("Protein Fusion: "):
            self.__isFusion = True
            s = self.__description
            self.__fusionGene = s[s.find("{")+1:s.find("}")]
            # check if fusion is defined in the correct format Gene1:Gene2
            try:
                self.__fusionPartner1, self.__fusionPartner2 = self.__fusionGene.split(":")
            except ValueError:
                raise IncorrectDescriptionFormat()

            # Check if fusion partners match the genes provided as inputs
            for gene in (self.__fusionPartner1, self.__fusionPartner2):
                if gene not in (self.__gene1, self.__gene2):
                    raise FusionGeneConflict(gene)
            if self.__fusionPartner1 == self.__gene1:
                self.__fusionTranscript1 = self.__transcript1
                self.__fusionTranscript2 = self.__transcript2
            else:
                self.__fusionTranscript2 = self.__transcript1
                self.__fusionTranscript1 = self.__transcript2

            # Is fusion known?
            if self.__fusionGene in oncoKb:
                self.__isKnownFusion = True
            else:
                self.__isKnownFusion = False

        # IF NOT FUSION
        else:
            self.__fusionGene = None
            self.__isFusion = False
            self.__isKnownFusion = False


        

    # define breakpoint order   (NOT USED AT THIS TIME) 
    def __set_bkpOrder(self):
        self.__bkpOrder = {}
        self.__tranlocationCoordOrder = {}

        # set annotation order based on fusion on not
        if self.__isFusion and (self.__isGene1InPanel or self.__isGene2InPanel):
            self.__bkpOrder[1]['gene'], self.__bkpOrder[2]['gene'] =\
                self.__fusionGene.split(":")
            self.__bkpOrder[1]['bkp'] = self.__gene_bkp_dict[self.__bkpOrder[1]['gene']]
            self.__bkpOrder[2]['bkp'] = self.__gene_bkp_dict[self.__bkpOrder[2]['gene']]
        elif self.__isGene1InPanel and self.__isGene2InPanel:
            self.__bkpOrder[1]['gene'], self.__bkpOrder[2]['gene'] =\
            self.__gene1, self.__gene2
            self.__bkpOrder[1]['bkp'], self.__bkpOrder[2]['bkp'] =\
                self.__bkp1, self.__bkp2
        elif self.__isGene1InPanel:
            self.__bkpOrder[1]['gene'] = self.__gene1
            self.__bkpOrder[1]['bkp'] = self.__bkp2
        elif self.__isGene1InPanel:
            self.__bkpOrder[2]['gene'] = self.__gene2
            self.__bkpOrder[2]['bkp'] = self.__bkp2
        else:
            raise GenesNotInPanel()
        
        # set translocation coordinate order
        if self.__svtype == "TRANSLOCATION":
        try:
            chr1, pos1 = self.__bkp1.split(":")
            chr2, pos2 = self.__bkp2.split(":")
        except ValueError:
            raise IncorrectBkpFormat()
        
        if chr1 == "X":
            self.__tranlocationCoordOrder[1] = self.__bkp1
            self.__tranlocationCoordOrder[2] = self.__bkp2
        elif chr2 == "X":
            self.__tranlocationCoordOrder[1] = self.__bkp2
            self.__tranlocationCoordOrder[2] = self.__bkp1  
        elif int(chr1) < int(chr2):
            self.__tranlocationCoordOrder[1] = self.__bkp1
            self.__tranlocationCoordOrder[2] = self.__bkp2
        elif int(chr1) > int(chr2):
            self.__tranlocationCoordOrder[1] = self.__bkp2
            self.__tranlocationCoordOrder[2] = self.__bkp1
        else:
            raise BreakpointConflictForTranslocation()
        
# not used at the moment
class gene(object):
    def __init__(self, geneID, bkp):
        self.id = geneID
        self.bkp = bkp

    def __expand(self, target_panel, panelKinase, refFlat):
        try:
            self.chr, self.pos = self.bkp.split(":")
        except ValueError:
            raise IncorrectBkpFormat()

        if refFlat[self.id]:
            self.transcript = refFlat[self.id]
        else:
            raise CanonicalTranscriptNotFound(self.id)

        if target_panel[self.id]:
            self.isPanel = True
        else:
            self.isPanel = False

        if panelKinase[self.id]:
            self.isKinase = True
        else:
            self.isKinase = False

            
