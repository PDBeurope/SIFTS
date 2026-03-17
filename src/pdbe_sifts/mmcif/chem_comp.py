#!/usr/bin/env python3

import pickle as pickle
from collections.abc import Mapping
from importlib.resources import files

from pdbe_sifts.base.log import logger

STANDARD_AA = {
    "ALA": "A",
    "ARG": "R",
    "ASP": "D",
    "ASN": "N",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


class ChemCompMapping:
    _dictionary: Mapping[str, str] = {}

    def __init__(self):
        cc_file = files("pdbe_sifts.data").joinpath(
            "three_to_one_letter_mapping.csv"
        )

        self.hydrate(cc_file)

    def hydrate(self, cc_file):
        logger.debug(f"CC_DICT: {cc_file}")

        logger.debug(f"Loading chem_comp from mapping file {cc_file}")
        with open(cc_file) as f:
            for line in f:
                three, one_letter = line.strip().split(",")
                self._dictionary[three] = one_letter

        if not self._dictionary:
            raise ValueError(
                f"Error in chemp_comp file/pkl, its empty!Please check file {cc_file}"
            )

        logger.debug(f"Loaded {len(self._dictionary)} mapping entries")
        self.__dict__ = self._dictionary

    def _check_chem_comp(self):
        """
        Sanity Check to ensure the chemp_comp always have correct value of 20 a.a
        """
        not_found = []
        for res in STANDARD_AA:
            if (
                res in self._dictionary
                and self._dictionary[res] == STANDARD_AA[res]
            ):
                continue
            else:
                if res in self._dictionary:
                    my_key = res
                    my_val = self._dictionary[res]
                else:
                    my_key = None
                    my_val = None

                not_found.append((res, STANDARD_AA[res], "vs", my_key, my_val))

        if not_found:
            logger.info("Not found:", not_found)
            raise Exception(
                f"""Error, no/mis-match of standard amino acid in your CC
                Please check {self.CC_DICT}
                """
            )

    def get(self, three):
        """Get one letter mapping from three letter."""
        if isinstance(three, list):
            out = ""
            for t in three:
                out += self.get(t)
            return out
        try:
            return self.__dict__.get(three.upper(), "X")
        except AttributeError:
            logger.debug(f"Cannot map for 3-letter code {three}")
            return
