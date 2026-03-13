class BatchRunException(Exception):
    """Raised when a batch run fails (e.g. no entries to process)."""

    pass


class ReleaseCheckFailedException(Exception):
    def __init__(self, entry_id, message, context=None):
        self.entry_id = entry_id
        self.context = context
        self.message = message
        super().__init__(message)


class ProcessFailedError(Exception):
    pass


class EntryFailedException(Exception):
    pass


class EntryTimedOutException(Exception):
    pass


class TooManyFailedEntries(Exception):
    """Thrown when entries that have failed are beyond an expected threshold"""

    pass


class ObsoleteUniProtError(Exception):
    """Exception raised when a Uniprot entry is obsolete."""

    pass


class AccessionNotFound(Exception):
    """Exception raised when a Uniprot entry does not exist.

    Uniprot API should return a 400 or 404 for this to be raised
    """

    pass


class SplitAccessionError(Exception):
    """Raised when an accession cannot be split or resolved during segment generation."""

    pass


class NoSegmentsError(Exception):
    """Raised when no SIFTS segments are written to mmCIF output."""

    pass


class NotAPolyPeptide(Exception):
    """Raised when a mmCIF entry contains no polypeptide chain."""

    pass
