#!/usr/bin/env python3
import itertools
from collections.abc import Iterable, Mapping
from operator import itemgetter

from pdbe_sifts.base.log import logger
from pdbe_sifts.segments_generation.alignment.helper import SMapping
from pdbe_sifts.unp.unp import UNP

# Parameterised query — the caller passes the entry value as a DuckDB `?` parameter
# to prevent SQL injection (CWE-89). Never interpolate entry directly into SQL.
_SELECTION_QUERY = """
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
        WHERE lower(entry) = ?
    ) sub
    WHERE sub.rn = 1
    ORDER BY entity
"""


def get_curated_db_mappings(
    pdbid: str,
    chains: Iterable,
    conn,
    chain_to_entity: Mapping[str, str],
) -> dict[str, list[SMapping]]:
    """Get mappings from database.

    Args:
        pdbid(str): PDB ID Code
        chains(Iterable): list of chains
        conn: Database connection

    Returns:
        Mapping[str, List[SMapping]]: Mapping for each chain of the sequence
    """
    rows = conn.execute(_SELECTION_QUERY, [pdbid])
    entity_mappings: Mapping[str, list[SMapping]] = {}
    for key, mapping in itertools.groupby(
        iter(rows.fetchall()), key=itemgetter(0, 1)
    ):
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
            logger.warning(
                f"Database returned a chain {chain} that's not a poly?"
            )
            del mappings[chain]

    return mappings
