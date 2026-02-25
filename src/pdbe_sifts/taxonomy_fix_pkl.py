#!/usr/bin/env python3

import os
import pickle as pickle
from collections.abc import Mapping
from importlib.resources import files

class TaxonomyFix:
    _dictionary: Mapping[str, Mapping[str, str]] = {}

    def __init__(self, conn=None, pkl_file=None):
        pkl_file = files("pdbe_sifts.data").joinpath("nf90_taxid.pkl")
        self._dictionary = {}

        if os.path.isfile(pkl_file):
            with open(pkl_file, "rb") as f:
                self._dictionary = pickle.load(f)

    def get(self, entry, entity):
        try:
            taxid = self.__dict__[entry][entity]
            return taxid
        except KeyError:
            return None

