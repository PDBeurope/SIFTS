#!/usr/bin/env python3
import itertools
from collections.abc import Iterable, Mapping
from operator import itemgetter

from pdbe_sifts.base.log import logger
# from pdbe_sifts.config import Config
from pdbe_sifts.segments_generation.alignment.helper import SMapping
from pdbe_sifts.unp.unp import UNP

# conf = Config()


def get_selection_queries(entry):
    sql = f"""
        SELECT
            entry, entity, accession, target_start, target_end
        FROM hits
        WHERE lower(entry) = '{entry}'
        and hit_rank = 1
        order by entry, entity, accession, target_start, target_end
    """
    return sql


def get_curated_db_mappings(pdbid, entities: Iterable, conn, unp_dir):
    """Get mappings from database.

    Args:
        pdbid(str): PDB ID Code
        entities(Iterable): list of entities
        conn: Database connection

    Returns:
        Mapping[str, List[SMapping]]: Mapping for each entity
    """
    sql = get_selection_queries(pdbid)

    rows = conn.execute(sql)
    mappings: Mapping[str, list[SMapping]] = {}
    for key, mapping in itertools.groupby(iter(rows.fetchall()), key=itemgetter(0, 1)):
        _, entity = key
        entity_mapping: list[SMapping] = []
        for _, _, accession, range_start, range_end in mapping:
            unp = UNP(accession, unp_dir)
            smap = SMapping(unp.accession, range_start, range_end)
            entity_mapping.append(smap)

        mappings[entity] = entity_mapping

    # chains_from_db = list(mappings.keys())
    # for chain in chains_from_db:
    #     if chain not in chains:
    #         logger.warning(f"Database returned a chain {chain} that's not a poly?")
    #         del mappings[chain]
    return mappings
