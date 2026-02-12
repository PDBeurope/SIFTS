#!/usr/bin/env python3

from pathlib import Path
from gemmi import cif
from pdbe_sifts.mmcif.chem_comp import ChemCompMapping


def extract_entities(cif_file):
    """
    Extract entity information from an mmCIF file.
    
    Args:
        cif_file: Path to the mmCIF file
        
    Returns:
        dict: Dictionary mapping entity_id to tuple (sequence, taxonomy_id)
              Example: {'1': ('MKLPQ...', 9606), '2': ('AATVR...', None)}
              
    Raises:
        FileNotFoundError: If the mmCIF file does not exist
    """
    if not Path(cif_file).exists():
        raise FileNotFoundError(f"The mmcif file {cif_file} does not exist.")
    
    block = cif.read(cif_file).sole_block()
    cc = ChemCompMapping()
    
    # Get required categories
    poly_seq = block.get_mmcif_category("_pdbx_poly_seq_scheme")
    entity_poly = block.get_mmcif_category("_entity_poly")
    src_nat = block.get_mmcif_category("_entity_src_nat")
    src_gen = block.get_mmcif_category("_entity_src_gen")
    src_syn = block.get_mmcif_category("_pdbx_entity_src_syn")
    
    if not poly_seq:
        return {}
    
    # Get all entity IDs from poly_seq
    entities = poly_seq["entity_id"]
    unique_entities = sorted(set(entities))
    
    result = {}
    
    for entity_id in unique_entities:
        # Check if it's a polypeptide
        if not _is_polypeptide(entity_poly, entity_id):
            continue
            
        # Build sequence
        sequence = _build_sequence(poly_seq, entity_id, cc)
        
        # Get taxonomy ID
        tax_id = _get_taxonomy_id(entity_id, src_nat, src_gen, src_syn)
        
        result[entity_id] = (sequence, tax_id)
    
    return result


def _is_polypeptide(entity_poly, entity_id):
    """Check if the entity is a polypeptide."""
    if not entity_poly:
        return False
        
    for idx, entity in enumerate(entity_poly["entity_id"]):
        if entity == entity_id:
            if isinstance(entity_poly["type"], list):
                return "peptide" in entity_poly["type"][idx].lower()
            else:
                return "peptide" in entity_poly["type"].lower()
    
    return False


def _build_sequence(poly_seq, entity_id, cc):
    """Build the sequence for a given entity."""
    entities = poly_seq["entity_id"]
    seq = {}
    
    for idx, e in enumerate(entities):
        if e != entity_id:
            continue
            
        seq_id = int(poly_seq["seq_id"][idx])
        
        if seq_id not in seq:
            three_letter = poly_seq["mon_id"][idx]
            seq[seq_id] = cc.get(three_letter)
    
    return "".join([seq[key] for key in sorted(seq.keys())])


def _get_taxonomy_id(entity_id, src_nat, src_gen, src_syn):
    """Get taxonomy ID for an entity from various sources."""
    # Try natural source first
    tax_id = _extract_tax_from_category(src_nat, "pdbx_ncbi_taxonomy_id", entity_id)
    if tax_id is not None:
        return tax_id
    
    # Try genetically manipulated source
    tax_id = _extract_tax_from_category(src_gen, "pdbx_gene_src_ncbi_taxonomy_id", entity_id)
    if tax_id is not None:
        return tax_id
    
    # Try synthetic source
    tax_id = _extract_tax_from_category(src_syn, "ncbi_taxonomy_id", entity_id)
    
    return tax_id


def _extract_tax_from_category(category, field_name, entity_id):
    """Extract taxonomy ID from a specific category."""
    if not category:
        return None
        
    entities = category["entity_id"]
    
    if not isinstance(entities, list):
        entities = [entities]
        is_list = False
    else:
        is_list = True
    
    for idx, ent in enumerate(entities):
        if entity_id == ent:
            if is_list:
                value = category[field_name][idx]
            else:
                value = category[field_name]
            
            # Handle comma-separated lists (take first value)
            if value and "," in value:
                value = value.split(",")[0]
            
            try:
                return int(value)
            except (ValueError, TypeError):
                return None
    
    return None