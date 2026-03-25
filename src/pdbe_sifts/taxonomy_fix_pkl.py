#!/usr/bin/env python3

import gzip
import os
import pickle as pickle
from collections.abc import Mapping
from importlib.resources import files


class TaxonomyFix:
    """Lookup table for correcting taxonomy IDs using a bundled pickle file.

    The correction data is loaded from ``pdbe_sifts/data/nf90_taxid.pkl.gz``
    at instantiation time.  Lookups are performed by PDB entry and entity
    identifier.

    Attributes:
        _dictionary (Mapping[str, Mapping[str, str]]): Nested mapping of
            ``{entry: {entity: taxid}}`` loaded from the pickle file.
    """

    _dictionary: Mapping[str, Mapping[str, str]] = {}

    def __init__(self, conn=None, pkl_file=None):
        """Load the taxonomy-fix dictionary from the bundled pickle file.

        Args:
            conn: Unused; reserved for future database-backed implementations.
            pkl_file: Unused; the path is always resolved from the package
                data directory (``pdbe_sifts/data/nf90_taxid.pkl.gz``).
        """
        pkl_file = files("pdbe_sifts.data").joinpath("nf90_taxid.pkl.gz")
        self._dictionary = {}

        if os.path.isfile(pkl_file):
            with gzip.open(pkl_file, "rb") as f:
                self._dictionary = pickle.load(f)

    def get(self, entry, entity):
        """Return the corrected taxonomy ID for a PDB entry/entity pair.

        Args:
            entry (str): The PDB entry identifier.
            entity (str | int): The entity identifier within the entry.

        Returns:
            str | None: The corrected taxonomy ID, or ``None`` if no
            correction is recorded for the given entry/entity combination.
        """
        try:
            taxid = self.__dict__[entry][entity]
            return taxid
        except KeyError:
            return None
