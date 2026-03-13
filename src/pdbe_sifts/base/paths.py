"""Path helpers for SIFTS: local data caches."""
from pathlib import Path

from pdbe_sifts.config import load_config

conf = load_config()


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
