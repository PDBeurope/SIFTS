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
        1. sifts_score          - composite score (adjusted + tax + dataset)
        2. pdb_cross_references - prefer accessions with more PDB structures
        3. adjusted_score       - alignment quality (identity x coverage)
        4. dataset_score        - Swiss-Prot > TrEMBL
        5. accession ASC        - deterministic alphabetic tiebreaker
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


def get_curated_db_mappings(pdbid, chains: Iterable, conn, chain_to_entity: Mapping[str, str]):
    """Get mappings from database.

    Args:
        pdbid(str): PDB ID Code
        chains(Iterable): list of chains
        conn: Database connection

    Returns:
        Mapping[str, List[SMapping]]: Mapping for each chain of the sequence
    """
    sql = get_selection_queries(pdbid)

    rows = conn.execute(sql)
    entity_mappings: Mapping[str, list[SMapping]] = {}
    for key, mapping in itertools.groupby(iter(rows.fetchall()), key=itemgetter(0, 1)):
        _, entity = key
        entity_mapping: list[SMapping] = []
        for _, _, accession, target_start, target_end in mapping:
            try:
                if not accession:
                    continue
                unp = UNP(accession)
            except Exception:
                continue
            try:
                smap = SMapping(unp.ad_dbref_auto_acc, target_start, target_end)
            except Exception:
                smap = SMapping(unp.accession, target_start, target_end)
            entity_mapping.append(smap)

        entity_mappings[entity] = entity_mapping

    # Build reverse mapping: entity_id -> [chains that share it]
    entity_to_chains: dict[str, list[str]] = {}
    for chain in chains:
        entity_id = str(chain_to_entity.get(chain, ""))
        if entity_id:
            entity_to_chains.setdefault(entity_id, []).append(chain)

    # Expand entity mappings to all chains sharing that entity
    mappings: dict[str, list[SMapping]] = {}
    for entity_id, smappings in entity_mappings.items():
        for chain in entity_to_chains.get(entity_id, []):
            mappings[chain] = smappings
        if entity_id not in entity_to_chains:
            logger.warning(
                f"DB returned entity {entity_id} for {pdbid} with no matching poly chain"
            )

    # Warn about DB chains not in our poly chain list
    for chain in list(mappings.keys()):
        if chain not in chains:
            logger.warning(f"Database returned a chain {chain} that's not a poly?")
            del mappings[chain]

    return mappings
