class Error(Exception):
    '''Base class for other exceptions'''
    pass


class MissingCytoBandError(Error):
    '''Raised when no cytoband was identified for a breakpoint'''

    def __init__(self):
        Exception.__init__(self, "No cytobands identified for the breakpoint.")


class MultipleCytoBandError(Error):
    '''Raised when multiple cytoband were identified for a breakpoint'''

    def __init__(self):
        Exception.__init__(
            self, "Multiple cytobands identified for the breakpoint.")


class CanonicalTranscriptNotFound(Error):
    '''Raised when canonical transcript not found'''

    def __init__(self, gene):
        Exception.__init__(
            self, "Cannot find a canonical transcript for the gene: " + gene)


class IncorrectGenesFormat(Error):
    '''Raised when gene1 / gene2 format does not match expected'''

    def __init__(self):
        Exception.__init__(
            self, "Input genes are not in required format. Should be in 'Gene1 / Gene2' format.")


class IncorrectBkpFormat(Error):
    '''Raised when a bkp format does not match expected'''

    def __init__(self, bkp):
        Exception.__init__(
            self, "Format of breakpoint " + str(bkp) + "does not match the required format. Example '7:1234567'.)


class IncorrectBkpFormat(Error):
    '''Raised when a bkp format does not match expected'''

    def __init__(self, bkp):
        Exception.__init__(
            self, "Format of breakpoint " + str(bkp) + "does not match the required format. Example '7:1234567'.)
#define reamining exception


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
            self, "Cannot determine cDNA annotation for the breakpoint: " + bkp
        )