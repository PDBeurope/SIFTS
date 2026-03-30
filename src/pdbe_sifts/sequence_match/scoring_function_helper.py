from ete4 import NCBITaxa

from pdbe_sifts.base.log import logger


def get_tax_weight(query_taxid: int, target_taxid: int) -> int:
    """
    Heuristic taxonomic weight based on the closest common ancestor (CCA).

    Rules:
      - same taxid -> 200
      - CCA distance from query: 0 -> 100, 1 -> 50, 2 -> 25, else -> 0

    Notes:
      - We compute the *deepest* shared ancestor (closest to the leaves), not the root.
    """
    if query_taxid == target_taxid:
        return 200

    ncbi = NCBITaxa()

    try:
        lineage1 = ncbi.get_lineage(query_taxid)[::-1]  # leaf → root
        lineage2 = ncbi.get_lineage(target_taxid)[::-1]  # leaf → root
    except Exception as e:
        logger.debug(f"Error computing lineage: {e}")
        return 0

    # First common ancestor in leaf→root order is the deepest (closest) shared node
    common_ancestors = [tax for tax in lineage1 if tax in lineage2]
    if not common_ancestors:
        return 0

    dist1 = lineage1.index(common_ancestors[0])
    return {0: 100, 1: 50, 2: 25}.get(dist1, 0)
