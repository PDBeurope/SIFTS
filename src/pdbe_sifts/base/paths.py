"""Path helpers for SIFTS: local data caches and configuration accessors.

`load_config()` is called once here. All other modules must import the getters
below instead of calling `load_config()` themselves.
"""
from pathlib import Path

from pdbe_sifts.config import load_config

conf = load_config()


# ── Cache path helpers ────────────────────────────────────────────────────────

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


# ── Configuration getters ─────────────────────────────────────────────────────

def get_conf_user_base_dir():
    return conf.user.base_dir


def get_conf_user_target_db():
    return conf.user.target_db


def get_conf_data_entry_dir():
    return conf.location.work.data_entry_dir


def get_conf_mmseqs_sensitivity() -> float:
    return conf.alignment.mmseqs.sensitivity


def get_conf_mmseqs_min_seq_id() -> float:
    return conf.alignment.mmseqs.min_seq_id


def get_conf_mmseqs_alignment_mode() -> int:
    return conf.alignment.mmseqs.alignment_mode


def get_conf_mmseqs_db_load_mode() -> int:
    return conf.alignment.mmseqs.db_load_mode


def get_conf_blastp_evalue() -> float:
    return conf.alignment.blastp.evalue


