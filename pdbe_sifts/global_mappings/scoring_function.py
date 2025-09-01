import argparse
import io
import csv
import json
from timeit import default_timer as timer
from pathlib import Path

from ete4 import NCBITaxa

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import get_date, make_path
from pdbe_sifts.unp.unp import UNP

# def base_score(identity, coverage):
#     return identity * coverage

def mismatch_penalty(nb_mismatches, seq_len):
    return 1.0 * nb_mismatches / seq_len

def adjusted_score(identity, coverage, nb_mismatches, seq_len):
    if coverage>0.7:
        return 1.5 * ((identity*coverage)*(1-mismatch_penalty(nb_mismatches, seq_len))) * 1000
    else:
        return (identity*coverage)*(1-mismatch_penalty(nb_mismatches, seq_len)) * 1000

def get_tax_weight(query_taxid: int, target_tax_id: int) -> int:
    if query_taxid == target_tax_id:
        return 200

    ncbi = NCBITaxa()
    
    try:
        lineage1 = ncbi.get_lineage(query_taxid)[::-1] # taxid --> root
        lineage2 = ncbi.get_lineage(target_tax_id)[::-1] # taxid --> root
    except ValueError as e:
        logger.warning(f'Taxonomy ID is not valid. Error: {e}')
        return 0

    common_ancestors = [tax for tax in lineage1 if tax in lineage2]

    for val in common_ancestors:
        dist1 = lineage1.index(val)
        if dist1==1:
            return 100
        elif dist1==2:
            return 50
        elif dist1==3:
            return 25
        else:
            return 0

def get_unp_info(accession, unp_dir):
    unp_obj = UNP(accession, unp_dir=unp_dir)
    dataset_score = 10 if unp_obj.dataset=='Swiss-Prot' else -10
    ref_prot_score = 100 if unp_obj.keywords and 'REFERENCE PROTEOME' in unp_obj.keywords else 0
    pdb_references_number_score = len(unp_obj.dbreferences.get('PDB', [])) * 0.1
    return {'dataset_score': dataset_score,
            'ref_prot_score': ref_prot_score, 
            'n_pdb_score': pdb_references_number_score}

def get_final_score(mapping, unp_dir):
    accession = mapping['accession']
    query_tax_id = mapping['query_tax_id']
    target_tax_id = mapping['target_tax_id']
    identity = mapping['identity']
    coverage = mapping['coverage']
    mismatch = mapping['mismatch']
    qlen = mapping['query_len']
    adj_score = adjusted_score(identity, coverage, mismatch, qlen)
    tax_weight = get_tax_weight(query_tax_id, int(target_tax_id))
    try:
        unp_data = get_unp_info(accession, unp_dir)
        dataset_score = unp_data['dataset_score']
        ref_prot_score = unp_data['ref_prot_score']
        n_pdb_score = unp_data['n_pdb_score']
    except ValueError as e:
        dataset_score = 0
        ref_prot_score = 0
        n_pdb_score = 0
        logger.warning(f'{e}')
    return adj_score + tax_weight + dataset_score + ref_prot_score + n_pdb_score


def get_ranked_mappings(mappings, unp_dir):
    ranked_mappings = {}
    for entry in mappings:
        ranked_mappings[entry] = {}
        for entity in mappings[entry]:
            ranked_mappings[entry][entity] = []
            scored_mappings = []
            for i, mapping in enumerate(mappings[entry][entity]):
                score = get_final_score(mapping, unp_dir)
                ind = i
                scored_mappings.append((i, score))
            sorted_list = sorted(scored_mappings, key=lambda x: x[1], reverse=True)
            for tpl in sorted_list:
                ind = tpl[0]
                score = tpl[1]
                mapping = mappings[entry][entity][ind]
                mapping['sifts_score'] = score
                ranked_mappings[entry][entity].append(mapping)
    return ranked_mappings
    
