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


class CanonicalTranscriptError(Error):
    '''Raised when canonical transcript not found'''

    def __init__(self):
        Exception.__init__(
            self, "Cannot find a canonical transcript for the breakpoint.")


class MultipleCanonicalTxError(Error):
    '''Raised when multiple canonical transcript versions exist'''

    def __init__(self):
        Exception.__init__(
            self, "Multiple versions of canonical transcript exists.")


#define reamining exception