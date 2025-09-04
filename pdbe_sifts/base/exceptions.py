class ObsoleteUniProtError(Exception):
    """Exception raised when a Uniprot entry is obsolete."""

    pass


class AccessionNotFound(Exception):
    """Exception raised when a Uniprot entry does not exist.

    Uniprot API should return a 400 or 404 for this to be raised
    """

    pass
