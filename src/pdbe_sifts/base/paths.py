"""Path helpers for SIFTS: local data caches and configuration accessors.

`load_config()` is called once here. All other modules must import the getters
below instead of calling `load_config()` themselves.
"""

from pathlib import Path

from pdbe_sifts.config import load_config

conf = load_config()

UNIPROT_PDB_TSV_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz"


def uniprot_cache_dir(accession: str, base_dir: str | None = None) -> str:
    """Local cache directory for a UniProt accession.

    Structure: {base}/{acc[0]}/{acc}/
    Defaults to conf.cache.uniprot when base_dir is not provided.
    """
    base = base_dir or conf.cache.uniprot
    return str(Path(base, accession[0], accession))


def ccd_cache_path(three_letter_code: str, base_dir: str | None = None) -> str:
    """Local CIF file path for a CCD chemical component.

    Structure: {base}/{code[0]}/{code}/{code}.cif
    Defaults to conf.cache.ccd when base_dir is not provided.
    """
    code = three_letter_code.upper()
    base = base_dir or conf.cache.ccd
    return str(Path(base, code[0], code, f"{code}.cif"))


def get_conf_user_base_dir() -> str | None:
    """Return the user-configured base output directory, or ``None`` if unset."""
    return conf.user.base_dir


def get_conf_user_target_db() -> str | None:
    """Return the user-configured target sequence database path, or ``None`` if unset."""
    return conf.user.target_db


def get_conf_data_entry_dir() -> str | None:
    """Return the configured per-entry data directory path, or ``None`` if unset."""
    return conf.location.work.data_entry_dir


def get_conf_mmseqs_sensitivity() -> float:
    """Return the MMseqs2 sensitivity parameter from configuration."""
    return conf.alignment.mmseqs.sensitivity


def get_conf_mmseqs_min_seq_id() -> float:
    """Return the MMseqs2 minimum sequence identity threshold from configuration."""
    return conf.alignment.mmseqs.min_seq_id


def get_conf_mmseqs_alignment_mode() -> int:
    """Return the MMseqs2 alignment mode integer from configuration."""
    return conf.alignment.mmseqs.alignment_mode


def get_conf_mmseqs_db_load_mode() -> int:
    """Return the MMseqs2 database load mode integer from configuration."""
    return conf.alignment.mmseqs.db_load_mode


def get_conf_mmseqs_prefilter_mode() -> int:
    """Return the MMseqs2 prefilter mode integer from configuration."""
    return conf.alignment.mmseqs.prefilter_mode


def get_conf_mmseqs_index_subset() -> int:
    """Return the MMseqs2 index subset integer from configuration."""
    return conf.alignment.mmseqs.index_subset


def get_conf_blastp_evalue() -> float:
    """Return the BLASTP e-value cutoff from configuration."""
    return conf.alignment.blastp.evalue


def get_conf_unp_pdb_xrefs_path() -> Path | None:
    """Return the path to the uniprot_pdb DuckDB index from config, or None."""
    val = conf.user.unp_pdb_xrefs
    return Path(val) if val else None
