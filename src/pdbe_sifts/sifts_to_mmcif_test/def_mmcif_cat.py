#!/usr/bin/env python3
"""
@following categories/items written in sifts.cif.gz
"""

SIFTS_NEW_CAT = [
    "atom_site",
    "pdbx_sifts_unp_segments",
    "pdbx_sifts_xref_db_segments",
    "pdbx_sifts_xref_db",
]
SIFTS_ATOMSITE_ITEM = [
    "id",
    "pdbx_label_index",
    "pdbx_sifts_xref_db_name",
    "pdbx_sifts_xref_db_acc",
    "pdbx_sifts_xref_db_num",
    "pdbx_sifts_xref_db_res",
]


# New mmcif categories introduced with their respective items

NEW_MMCIF_CAT = {
    "_pdbx_sifts_unp_segments": [
        "entity_id",
        "asym_id",
        "unp_acc",
        "segment_id",
        "instance_id",
        "unp_start",
        "unp_end",
        "seq_id_start",
        "seq_id_end",
        "best_mapping",
        "identity",
    ],
    "_pdbx_sifts_xref_db_segments": [
        "entity_id",
        "asym_id",
        "xref_db",
        "xref_db_acc",
        "domain_name",
        "segment_id",
        "instance_id",
        "seq_id_start",
        "seq_id_end",
    ],
    "_pdbx_sifts_xref_db": [
        "entity_id",
        "asym_id",
        "seq_id_ordinal",
        "seq_id",
        "mon_id",
        "mon_id_one_letter_code",
        "unp_res",
        "unp_num",
        "unp_acc",
        "unp_segment_id",
        "unp_instance_id",
        "res_type",
        "observed",
        "mh_id",
        "xref_db_name",
        "xref_db_acc",
        "xref_domain_name",
        "xref_db_segment_id",
        "xref_db_instance_id",
    ],
}

# Primary keys info used to checking the output,
# primary keys should always be present and should not be null

PRIMARY_KEYS = {
    "_pdbx_sifts_unp_segments": [
        "asym_id",
        "entity_id",
        "instance_id",
        "segment_id",
        "unp_acc",
    ],
    "_pdbx_sifts_xref_db_segments": [
        "asym_id",
        "entity_id",
        "instance_id",
        "segment_id",
        "xref_db",
        "xref_db_acc",
    ],
    "_pdbx_sifts_xref_db": ["asym_id", "entity_id", "seq_id", "seq_id_ordinal"],
    "_atom_site": ["id"],
}


# Categotries to compare for tracking changes in mapping

DELTA_CAT = {
    "_pdbx_sifts_unp_segments": [
        "asym_id",
        "unp_acc",
        "unp_start",
        "unp_end",
        "seq_id_start",
        "seq_id_end",
    ],
    "_pdbx_sifts_xref_db_segments": [
        "asym_id",
        "xref_db",
        "xref_db_acc",
        "seq_id_start",
        "seq_id_end",
    ],
}

# Databases for which the mappings is written in mmcif

MAPPED_DB = ["UNP", "Pfam", "CATH", "SCOP2"]
