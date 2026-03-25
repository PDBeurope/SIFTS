#!/usr/bin/env python3


import pickle as pickle
from collections.abc import Mapping

from pdbe_sifts.base.log import logger
from pdbe_sifts.unp.unp import UNP

from . import mmcif_helper
from .chain import Chain


class Entry:
    """Top-level container for a PDB entry used by the SIFTS pipeline.

    Parses the mmCIF file, filters to polypeptide chains, and constructs
    :class:`~pdbe_sifts.mmcif.chain.Chain` objects for each one.

    Attributes:
        pdbid (str): Four-letter PDB identifier.
        mmcif (mmcif_helper.mmCIF): Parsed mmCIF data accessor.
        chains (Mapping[str, Chain]): Author chain ID → Chain object for
            every polypeptide chain in the entry.
        sequences (dict): Entity ID → one-letter sequence string.
        accessions (Mapping[str, UNP]): UniProt accession → UNP object
            cache, populated during alignment.
    """

    def __init__(self, pdbid, chem_comp_dict, cif_file: str):
        """Initialise an Entry by parsing the given mmCIF file.

        Reads all mmCIF categories required by the SIFTS pipeline, then
        constructs a :class:`Chain` object for every polypeptide chain
        found in ``_pdbx_poly_seq_scheme``.  Non-polypeptide chains are
        silently skipped.

        Args:
            pdbid: Four-letter PDB identifier (e.g. ``"1abc"``).
            chem_comp_dict: A :class:`~pdbe_sifts.mmcif.chem_comp.ChemCompMapping`
                instance used to convert three-letter residue codes to
                one-letter codes.
            cif_file: Absolute path to the ``.cif`` or ``.cif.gz`` file.

        Raises:
            FileNotFoundError: If *cif_file* does not exist on disk.
            NotAPolyPeptide: If the file contains no polypeptide sequence
                data (``_pdbx_poly_seq_scheme`` is absent).
        """
        self.pdbid = pdbid
        self.mmcif = mmcif_helper.mmCIF(pdbid, chem_comp_dict, cif_file)
        self.chains: Mapping[str, Chain] = {}

        # map entity -> sequence
        self.sequences = {}

        # map uniprot -> UNP object
        self.accessions: Mapping[str, UNP] = {}

        for key, val in list(self.mmcif.get_chains().items()):
            entity_id = val[0]
            struct_asym = val[1]

            # Reject non polypeptides
            if not self.mmcif.is_poly(entity_id):
                logger.debug(f"Rejecting chain {key}. Not a polypeptide")
                continue

            if entity_id not in self.sequences:
                self.sequences[val[0]] = self.mmcif.get_sequence(val[0])
            logger.debug(f"Adding {key} to chains")
            self.chains[key] = Chain(
                self.mmcif,
                self.pdbid,
                key,
                entity_id,
                struct_asym,
                self.sequences[entity_id],
                self,
            )
