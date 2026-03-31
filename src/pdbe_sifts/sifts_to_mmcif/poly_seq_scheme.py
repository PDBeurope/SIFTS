"""Backward-compatibility re-export — module has moved to ``pdbe_sifts.mmcif.poly_seq_scheme``."""

from pdbe_sifts.mmcif.poly_seq_scheme import (  # noqa: F401
    SCHEME_COLUMNS,
    build_pdbx_poly_seq_scheme,
    build_pdbx_poly_seq_scheme_all,
    get_all_entity_chain_pairs,
    get_canonical_residues,
    get_coordinate_residues,
    write_pdbx_poly_seq_scheme,
)
