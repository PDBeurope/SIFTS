from ete4 import NCBITaxa

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