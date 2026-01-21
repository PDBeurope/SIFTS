#!/usr/bin/env python3

import pickle as pickle
from collections.abc import Mapping

from pdbe_sifts.base.log import logger
from pdbe_sifts.config import Config
from pdbe_sifts.unp.unp import UNP

from . import mmcif_helper
from .chain import Chain
from .entity import Entity


class Entry:
    """Docstring for Entry."""

    def __init__(self, name, cif_file):
        # self.cif_dir = cif_dir or conf.location.work.data_entry_dir
        self.cif_file = cif_file

        self.name = name
        self.mmcif = mmcif_helper.mmCIF(self.name, self.cif_file)
        self.chains: Mapping[str, Chain] = {}
        self.entities: Mapping[str, Entity] = {}

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
                self.name,
                key,
                entity_id,
                struct_asym,
                self.sequences[entity_id],
                self,
            )
            if entity_id not in self.entities:
                self.entities[entity_id] = Entity(
                    self.mmcif,
                    self.name,
                    entity_id,
                    self,
                )
            self.entities[entity_id].add_chain(key, self.chains[key])
            

    def get_entity_seq_tax(self):
        entity_chain_seq = {}
        for key, chain in self.chains.items():
            ent = chain.entity_id
            tax = chain.tax_id
            seq = chain.sequence
            if ent not in entity_chain_seq:
                entity_chain_seq[ent] = (seq, tax)
        return entity_chain_seq