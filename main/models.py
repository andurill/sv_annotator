class bkp(object):
    def __init__(self, chrom, pos, gene, desc):
        try:
            self.chr = str(chrom)
            self.pos = int(pos)
            self.gene = str(gene)
            self.desc = str(desc)
        except TypeError:
            print("Could not create a new instance of bkp class due to inappropriate values for parameters.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

    def __expand(self, target_panel, panelKinase, refFlat):
        if refFlat[self.gene]:
            self.transcript = refFlat[self.gene]
        else:
            raise CanonicalTranscriptNotFound(self.gene)

        if target_panel[self.gene]:
            self.isPanel = True
        else:
            self.isPanel = False

        if panelKinase[self.gene]:
            self.isKinase = True
        else:
            self.isKinase = False


class sv(object):
    def __init__(self, svtype, bkp1, bkp2, genes, site1, site2, description, connection):
        self.__svtype = svtype
        self.__site1 = site1
        self.__site2 = site2
        self.__description = description
        self.__connection = connection
        try:
            self.__chr1, self.__pos1 = self.__bkp1.split(":") #IncorrectBkpFormat()
            self.__chr2, self.__chr2 = self.__bkp2.split(":") #IncorrectBkpFormat()
            self.__gene1, self.__gene2 = genes.split(" / ") #IncorrectGenesFormat()
        except ValueError:
            print("Could not create a new instance of sv class due to inappropriate values for parameters.")

    def __expand(self):
        self.bkp1 = bkp(self.__chr1, self.__pos1, self.__gene1, self.__site1)
        self.bkp2 = bkp(self.__chr2, self.__pos2, self.__gene2, self.__site2)

        # Define key variables
        if self.bkp1.InPanel is False and self.bkp2.InPanel is False:
            raise GenesNotInPanel()

        # Is SV intragenic?
        if self.bkp1.gene == self.bkp2.gene:
            self.__isIntragenic = True
        else:
            self.__isIntragenic = False

        # Fusion variables
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
            if self.__fusionPartner1 == self.bkp1.gene and self.__fusionPartner2 == self.bkp2.gene:
                self.__fusionPartner1, self.__fusionPartner2 = self.bkp1, self.bkp2
            elif self.__fusionPartner2 == self.bkp1.gene and self.__fusionPartner1 == self.bkp2.gene:
                self.__fusionPartner2, self.__fusionPartner1 = self.bkp1, self.bkp2      
            else:
                raise FusionGeneConflict()         

            # Is fusion known?
            if self.__fusionGene in oncoKb:
                self.__isKnownFusion = True
            else:
                self.__isKnownFusion = False
        else:
            self.__fusionGene = None
            self.__isFusion = False
            self.__isKnownFusion = False

        # Annotation variables
        if self.__svtype == "TRANSLOCATION":
            if self.bkp1.chrom.isdigit() and self.bkp2.chrom.isdigit():
                if int(self.bkp1.chrom) < int(self.bkp1.chrom):
                    self.__annotationPartner1, self.__annotationPartner2 = self.bkp1, self.bkp2
                else:
                    self.__annotationPartner1, self.__annotationPartner2 = self.bkp2, self.bkp1
            elif str(self.bkp1.chrom) == "X":
                self.__annotationPartner1, self.__annotationPartner2 = self.bkp1, self.bkp2
            elif str(self.bkp2.chrom) == "X":
                self.__annotationPartner1, self.__annotationPartner2 = self.bkp2, self.bkp1
            elif str(self.bkp1.chrom) == "Y":
                self.__annotationPartner1, self.__annotationPartner2 = self.bkp1, self.bkp2
            else:
                self.__annotationPartner1, self.__annotationPartner2 = self.bkp2, self.bkp1
        elif self.__isFusion is True:
            self.__annotationPartner1, self.__annotationPartner2 = self.__fusionPartner1, self.__fusionPartner2
        elif self.bkp1.isPanel is True and self.bkp2.isPanel is True:
            self.__annotationPartner1, self.__annotationPartner2 = self.bkp1, self.bkp2
        elif


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
