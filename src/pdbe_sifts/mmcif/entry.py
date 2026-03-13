#!/usr/bin/env python3


import pickle as pickle
from collections.abc import Mapping

from pdbe_sifts.base.log import logger
from pdbe_sifts.config import load_config
from pdbe_sifts.unp.unp import UNP

from . import mmcif_helper
from .chain import Chain

conf = load_config()


class Entry:
    """Docstring for Entry."""

    def __init__(self, pdbid, chem_comp_dict, cif_file: str):
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
