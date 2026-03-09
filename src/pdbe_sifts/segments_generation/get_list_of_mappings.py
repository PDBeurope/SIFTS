#!/usr/bin/env python3
import itertools
from collections.abc import Iterable, Mapping
from operator import itemgetter

from pdbe_sifts.base.log import logger
from pdbe_sifts.segments_generation.alignment.helper import SMapping
from pdbe_sifts.unp.unp import UNP


def get_selection_queries(entry):
    """Return the SQL query that fetches exactly one best hit per entity.

    Uses ROW_NUMBER() instead of DENSE_RANK (hit_rank) so that ties are
    broken deterministically and only a single accession is returned per
    entity — regardless of how many accessions share the top sifts_score.

    Tie-breaking order (descending quality):
        1. sifts_score          – composite score (adjusted + tax + dataset)
        2. pdb_cross_references – prefer accessions with more PDB structures
        3. adjusted_score       – alignment quality (identity × coverage)
        4. dataset_score        – Swiss-Prot > TrEMBL
        5. accession ASC        – deterministic alphabetic tiebreaker
    """
    sql = f"""
        SELECT entry, entity, accession, target_start, target_end
        FROM (
            SELECT
                entry, entity, accession, target_start, target_end,
                ROW_NUMBER() OVER (
                    PARTITION BY entry, entity
                    ORDER BY
                        sifts_score          DESC,
                        pdb_cross_references DESC,
                        adjusted_score       DESC,
                        dataset_score        DESC,
                        accession            ASC
                ) AS rn
            FROM hits
            WHERE lower(entry) = '{entry}'
        ) sub
        WHERE sub.rn = 1
        ORDER BY entity
    """
    return sql


def get_curated_db_mappings(pdbid, entities: Iterable, conn, unp_dir):
    """Get mappings from database.

    Args:
        pdbid(str): PDB ID Code
        entities(Iterable): list of entities
        conn: Database connection
        unp_dir: Directory for cached UniProt XML files

    Returns:
        Mapping[str, List[SMapping]]: Mapping for each entity
    """
    sql = get_selection_queries(pdbid)

    rows = conn.execute(sql)
    mappings: Mapping[str, list[SMapping]] = {}
    for key, mapping in itertools.groupby(iter(rows.fetchall()), key=itemgetter(0, 1)):
        _, entity = key
        entity_mapping: list[SMapping] = []
        for _, _, accession, target_start, target_end in mapping:
            unp = UNP(accession, unp_dir)
            smap = SMapping(unp.accession, target_start, target_end)
            entity_mapping.append(smap)
        mappings[entity] = entity_mapping

    return mappings
