import json
from typing import Dict, List, Tuple, Any

from ete4 import NCBITaxa
from pdbe_sifts.base.log import logger
from pdbe_sifts.unp.unp import UNP


def mismatch_penalty(nb_mismatches, seq_len):
    """Return a simple mismatch penalty proportional to the mismatch ratio."""
    if seq_len <= 0:
        return 0.0
    return 1.0 * nb_mismatches / seq_len

def adjusted_score(identity, coverage, nb_mismatches, seq_len):
    """
    Compute an adjusted alignment score combining identity, coverage and a mismatch penalty.
    """
    # if coverage>0.7:
    #     return 1.5 * ((identity*coverage)*(1-mismatch_penalty(nb_mismatches, seq_len))) * 1000
    # else:
    return (identity*coverage)*(1-mismatch_penalty(nb_mismatches, seq_len)) * 1000

def get_tax_weight(query_taxid: int, target_tax_id: int) -> int:
    """
    Heuristic taxonomic weight based on the closest common ancestor (CCA).

    Rules:
      - same taxid -> 200
      - CCA distance from query: 0 -> 100, 1 -> 50, 2 -> 25, else -> 0

    Notes:
      - We compute the *deepest* shared ancestor (closest to the leaves), not the root.
    """
    if query_taxid == target_tax_id:
        return 200

    ncbi = NCBITaxa()
    
    try:
        lineage1 = ncbi.get_lineage(query_taxid)[::-1] # taxid --> root
        lineage2 = ncbi.get_lineage(target_tax_id)[::-1] # taxid --> root
    except Exception as e:
        logger.debug(f'Error: {e}')
        return 0

    common_ancestors = [tax for tax in lineage1 if tax in lineage2]

    # Find deepest common ancestor by scanning lineage1 from leaf to root
    for val in common_ancestors:
        dist1 = lineage1.index(val)
        if dist1==0:
            return 100
        elif dist1==1:
            return 50
        elif dist1==2:
            return 25
        else:
            return 0

def get_unp_info(accession, unp_dir):
    """
    Retrieve UNP-derived scores for an accession.

    Returns:
        {
            "dataset_score": float,
            "n_pdb_score": int
        }
    """
    try:
        unp_obj = UNP(accession, unp_dir=unp_dir)
    except Exception as e:
        logger.debug(f'UNP load failed for {accession}:{e}')
        return {"dataset_score": 0.0, "pdb_references_number": 0}

    if unp_obj.dataset == "Swiss-Prot":
        if unp_obj.annotation_score == 0:
            dataset_score = 10
        else:
            dataset_score = 10 * unp_obj.annotation_score
    else:
        dataset_score = -50 + 10 * unp_obj.annotation_score

    pdb_references_number = len(unp_obj.dbreferences.get('PDB', []))
    return {'dataset_score': dataset_score,
            'pdb_references_number': pdb_references_number}

def get_final_score(mapping: dict, unp_dir: str) -> float:
    """
    Calculate the SIFTS score by combining multiple metrics:
      - Adjusted alignment score (identity, coverage, mismatches, length)
      - Taxonomic weight (query_tax_id vs. target_tax_id)
      - UNP-derived dataset score
    Also returns the number of PDB cross-references (not included in the score).

    Args:
        mapping: Dict with keys:
            'accession', 'query_tax_id', 'target_tax_id',
            'identity', 'coverage', 'mismatch', 'query_len'
        unp_dir: Path to the UNP data directory.

    Returns:
        (sifts_score: float, n_pdb: int)
    """
   # Extract data from mapping
    accession = mapping['accession']
    query_tax_id = mapping['query_tax_id']
    target_tax_id = mapping['target_tax_id']
    identity = mapping['identity']
    coverage = mapping['coverage']
    mismatch = mapping['mismatch']
    qlen = mapping['query_len']

    # Calculate base scores
    adj_score = adjusted_score(identity, coverage, mismatch, qlen)
    tax_weight = get_tax_weight(query_tax_id, int(target_tax_id))

    # Initialize default UNP scores
    unp_scores = {
        'dataset_score': 0,
        'pdb_references_number': 0
    }
    # Retrieve UNP scores if available
    try:
        unp_data = get_unp_info(accession, unp_dir)
        unp_scores.update({
            'dataset_score': unp_data.get('dataset_score', 0),
            'pdb_references_number': unp_data.get('pdb_references_number', 0),
        })
    except Exception as e:
        logger.debug(f"Failed to retrieve UNP data for {accession}: {e}")
    # Calculate sifts score
    sifts_score = (
        adj_score +
        tax_weight +
        unp_scores['dataset_score']
    )
    pdb_references_number = unp_scores['pdb_references_number']
    return sifts_score, pdb_references_number


def get_ranked_mappings(mappings, unp_dir):
    """
    Rank mappings by SIFTS score, per entry and per entity.

    Args:
        mappings: Nested dict structure:
            {
              entry: {
                entity: [ mapping_dict, ... ]
              },
              ...
            }
        unp_dir: Path to the UNP data directory.

    Returns:
        Same structure as input, but each mapping enriched with:
          - 'sifts_score': float
          - 'pdb_cross_references': int
        and lists sorted by 'sifts_score' descending.
    """
    ranked_mappings = {}
    for entry in mappings:
        ranked_mappings[entry] = {}
        for entity in mappings[entry]:
            ranked_mappings[entry][entity] = []
            scored_mappings = []
            for i, mapping in enumerate(mappings[entry][entity]):
                score, pdb_references_number = get_final_score(mapping, unp_dir)
                ind = i
                scored_mappings.append((i, score, pdb_references_number))
            sorted_list = sorted(scored_mappings, key=lambda x: x[1], reverse=True)
            for tpl in sorted_list:
                ind = tpl[0]
                score = tpl[1]
                pdb_references_number = tpl[2]
                mapping = mappings[entry][entity][ind]
                mapping['sifts_score'] = score
                mapping['pdb_cross_references'] = pdb_references_number
                ranked_mappings[entry][entity].append(mapping)
    return ranked_mappings
    
