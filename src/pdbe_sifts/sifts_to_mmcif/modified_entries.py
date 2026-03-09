"""
This script generates a list of entries that have been modified in SIFTS
for downstream processes in the release pipeline.

Also adds all new and modified entries this week to the list.
"""

from pdbe_sifts.base.log import logger
from pdbe_sifts.config import load_config

conf = load_config()


def run():
    entries_this_week = set()
    with open(conf.lists.entries_this_week) as f:
        for line in f:
            entries_this_week.add(line.strip())

    sifts_mapping_changes = set()
    with open(conf.lists.sifts_mapping_changes) as f:
        for line in f:
            sifts_mapping_changes.add(line.strip().split()[0])

    sifts_modified_entries = sorted(entries_this_week.union(sifts_mapping_changes))
    with open(conf.lists.sifts_modified_entries, "w+") as w:
        for entry in sifts_modified_entries:
            w.write(f"{entry}\n")

    logger.info(
        "{} entries modified this week written to {}".format(
            len(sifts_modified_entries), conf.lists.sifts_modified_entries
        )
    )
