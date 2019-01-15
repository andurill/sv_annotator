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
                self.__fusionPartner1, self.__fusionPartner2 = self.__fusionGene.split(
                    ":")
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



class Error(Exception):
    '''Base class for other exceptions'''
    pass


class MissingCytoBand(Error):
    '''Raised when no cytoband was identified for a breakpoint'''

    def __init__(self, bkp):
        Exception.__init__(self, "No cytobands identified for the breakpoint: " + bkp
        )


class MultipleCytoBand(Error):
    '''Raised when multiple cytoband were identified for a breakpoint'''

    def __init__(self):
        Exception.__init__(
            self, "Multiple cytobands identified for the breakpoint: " + bkp
            )


class CanonicalTranscriptNotFound(Error):
    '''Raised when canonical transcript not found'''

    def __init__(self, gene):
        Exception.__init__(
            self, "Cannot find a canonical transcript for the gene: " + gene
            )


class IncorrectGenesFormat(Error):
    '''Raised when gene1 / gene2 format does not match expected'''

    def __init__(self):
        Exception.__init__(
            self, "Input genes are not in required format. Should be in 'Gene1 / Gene2' format."
            )


class IncorrectBkpFormat(Error):
    '''Raised when a bkp format does not match expected'''

    def __init__(self, bkp):
        Exception.__init__(
            self, "Format of breakpoint " + str(bkp) + "does not match the required format. Example '7:1234567'".
            )


class IncorrectDescriptionFormat(Error):
    '''Raised when the input fusion description does not match expected'''

    def __init__(self):
        Exception.__init__(
            self, "Format of fusion description does not meet the required format. Example 'Protein Fusion: {FGFR3:TACC3}'"
        )


class GenesNotInPanel(Error):
    '''Raised when both genes not in panel'''

    def __init__(self):
        Exception.__init__(
            self, "Neither of the gene partners is in the target panel."
        )


class FusionGeneConflict(Error):
    '''Raised when one of the fusion parter does not match either of the input genes'''

    def __init__(self, gene):
        Exception.__init__(
            self, "The fusion partner, " + str(gene) + " does not match either of the two input genes."
        )


class cdnaNotFound(Error):
    '''Raised when vep cannot generate cdna annotation for a breakpoint'''

    def __init__(self, bkp):
        Exception.__init__(
            self, "Cannot determine cDNA annotation using vep for the breakpoint: " + bkp
        )


class IncorrectSVInputArguments(Error):
    '''Raised when the input arguments to sv class in incorrect'''

    def __init__(self):
        Exception.__init__(
            self, "The number arguments provided for sv class in correct. It should be a \
            comma-separater string of 'svtype, bkp1, bkp2, genes, site1, site2, description, connection'."
        )