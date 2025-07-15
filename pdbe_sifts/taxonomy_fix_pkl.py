#!/usr/bin/env python3

import os
import pickle as pickle
from collections.abc import Mapping

# from release.config import Config

# conf = Config()


class TaxonomyFix:
    _dictionary: Mapping[str, Mapping[str, str]] = {}

    def __init__(self, conn=None, pkl_file=None):
        conf = '/hps/software/users/pdbe/user/adamb/opensifts/opensifts/input_files/nf90_taxid.pkl'#pkl_file or conf.sifts.tax_fix_pkl

        if os.path.isfile(conf):
            with open(conf, "rb") as f:
                self._dictionary = pickle.load(f)
        elif conn:
            cursor = conn.cursor()
            cursor.execute("SELECT entry_id, entity_id, tax_id FROM taxonomy_fix")

            if cursor:
                for row in cursor:
                    entry_id, entity_id, tax_id = row[0], row[1], row[2]
                    self._dictionary.setdefault(entry_id, {})[entity_id] = tax_id

            cursor.close()

            with open(conf, "wb") as f:
                pickle.dump(self._dictionary, f, pickle.HIGHEST_PROTOCOL)
        self.__dict__ = self._dictionary

    def get(self, entry, entity):
        try:
            taxid = self.__dict__[entry][entity]
            return taxid
        except KeyError:
            return None

