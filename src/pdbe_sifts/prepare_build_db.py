from __future__ import annotations

import gzip
import re
import shutil
from pathlib import Path
from typing import TextIO

import requests

from pdbe_sifts.base.log import logger

DEFAULT_UNIPROT_SWISSPROT_FASTA_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/uniprot_sprot.fasta.gz"
)

_UNIPROT_FASTA_HEADER_RE = re.compile(
    r"^(?:sp|tr)\|(?P<accession>[^|]+)\|.*\bOX=(?P<tax_id>\d+)\b"
)


def prepare_build_db(
    output_fasta: str | Path,
    output_tax_mapping: str | Path,
    input_fasta: str | Path | None = None,
    force: bool = False,
) -> tuple[Path, Path]:
    """Prepare FASTA and taxonomy mapping inputs for target DB creation.

    When *input_fasta* is not provided, the reviewed UniProt Swiss-Prot FASTA
    is downloaded. The taxonomy mapping is parsed from UniProtKB FASTA headers
    with the shape ``>sp|ACCESSION|... OX=TAXID`` or
    ``>tr|ACCESSION|... OX=TAXID``.
    """
    fasta_path = _prepare_fasta(
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        force=force,
    )
    tax_mapping_path = Path(output_tax_mapping)
    write_uniprot_tax_mapping(fasta_path, tax_mapping_path)
    return fasta_path, tax_mapping_path


def _prepare_fasta(
    input_fasta: str | Path | None,
    output_fasta: str | Path,
    force: bool = False,
) -> Path:
    output_path = Path(output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if input_fasta:
        input_path = Path(input_fasta)
        if not input_path.exists():
            raise FileNotFoundError(f"Input FASTA not found: {input_path}")
        if input_path.resolve() == output_path.resolve():
            return output_path
        if output_path.exists() and not force:
            logger.info(
                "FASTA already present at %s; skipping copy", output_path
            )
            return output_path
        _copy_fasta(input_path, output_path)
        logger.info("Copied FASTA from %s to %s", input_path, output_path)
        return output_path

    if output_path.exists() and not force:
        logger.info(
            "FASTA already present at %s; skipping download", output_path
        )
        return output_path

    _download_file(DEFAULT_UNIPROT_SWISSPROT_FASTA_URL, output_path)
    return output_path


def write_uniprot_tax_mapping(
    fasta_path: str | Path,
    output_tax_mapping: str | Path,
) -> Path:
    """Write ``ACCESSION<TAB>TAXID`` rows parsed from a UniProtKB FASTA."""
    fasta = Path(fasta_path)
    output_path = Path(output_tax_mapping)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    skipped = 0
    with _open_text(fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            parsed = parse_uniprot_tax_header(line[1:].strip())
            if parsed is None:
                skipped += 1
                continue
            rows.append(parsed)

    if not rows:
        raise ValueError(
            "No UniProtKB taxonomy IDs were parsed from FASTA headers. "
            "Expected headers like '>sp|P12345|... OX=9606 ...' or "
            "'>tr|A0A000|... OX=9606 ...'."
        )

    with output_path.open("w", encoding="utf-8") as out_handle:
        for accession, tax_id in rows:
            out_handle.write(f"{accession}\t{tax_id}\n")

    logger.info(
        "Wrote %d taxonomy mapping rows to %s (%d headers skipped)",
        len(rows),
        output_path,
        skipped,
    )
    return output_path


def parse_uniprot_tax_header(header: str) -> tuple[str, str] | None:
    """Parse a UniProtKB FASTA header into ``(accession, tax_id)``."""
    match = _UNIPROT_FASTA_HEADER_RE.match(header)
    if match is None:
        return None
    return match.group("accession"), match.group("tax_id")


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _copy_fasta(input_path: Path, output_path: Path) -> None:
    input_is_gzip = input_path.name.endswith(".gz")
    output_is_gzip = output_path.name.endswith(".gz")

    if input_is_gzip == output_is_gzip:
        shutil.copyfile(input_path, output_path)
        return

    if output_is_gzip:
        with (
            _open_text(input_path) as input_handle,
            gzip.open(output_path, "wt", encoding="utf-8") as output_handle,
        ):
            shutil.copyfileobj(input_handle, output_handle)
        return

    with (
        _open_text(input_path) as input_handle,
        output_path.open("w", encoding="utf-8") as output_handle,
    ):
        shutil.copyfileobj(input_handle, output_handle)


def _download_file(url: str, output_path: Path) -> None:
    logger.info("Downloading %s to %s", url, output_path)
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        with output_path.open("wb") as out_handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    out_handle.write(chunk)
