#!/usr/bin/env python3

import requests
import shutil
from pathlib import Path
from datetime import datetime
from typing import Optional, List
from xml.etree import ElementTree
import funcy

from pdbe_sifts.base.log import logger

UNIPROT_REGEX = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
UNIPROT_API_BASE_URL = "https://www.uniprot.org/uniprot"

def get_uniprot_cache_dir(unp_accession, base_dir):
    """Path to UNP cache directory for specified UNP accession."""
    return str(Path(base_dir, unp_accession[0], unp_accession))


def fetch_uniprot_file(uniprot_id: str, file_type: str, unp_dir, fail_silently=False) -> str:
    """Fetches Uniprot file for a given Uniprot ID.

    First checks cache directory for the file. If not found, fetches from Uniprot.
    Cache is derived from config variable `cache.uniprot`.

    Args:
        uniprot_id (str): Uniprot ID to fetch.
        file_type (str): Type of file to fetch. One of: xml, json, fasta.
        fail_silently(bool): Return None on fail instead of raising error

    Raises:
        ValueError: If file_type is not one of xml, json, fasta.
        ObsoleteUniProtError: If the Uniprot entry is deleted (including Blank XML/FASTA).
        AccessionNotFound: If the Uniprot entry is not found (404 from Uniprot API)
    """
    try:

        if file_type not in ["xml", "json", "fasta"]:
            raise ValueError(
                f"Invalid file type {file_type}. Must be one of xml, json, fasta"
            )
        filename = f"{uniprot_id}.{file_type}"
        unp_file = Path(unp_dir, filename)

        unp_file.parent.mkdir(parents=True, exist_ok=True)

        if Path.exists(unp_file):
            logger.info(f"Fetched from cache: {unp_file}")
            return unp_file

        _unp_file_checks(file_type, filename, unp_file)

        return unp_file
    except Exception:
        if fail_silently:
            return
        raise


@funcy.retry(
    3,
    timeout=lambda a: 0.1**a,
    errors=(
        requests.exceptions.RequestException,
        requests.exceptions.ReadTimeout,
    ),
)

def _unp_file_checks(file_type, filename, unp_file):
    url = f"{UNIPROT_API_BASE_URL}/{filename}"
    logger.info(f"Fetching {url}")
    r = requests.get(url, timeout=5)

    if r.status_code in [404, 400]:
        raise AccessionNotFound(f"Uniprot entry {url} not found")

    r.raise_for_status()

    if file_type == "json":
        _check_json_response(r)
    elif file_type == "xml":
        _check_xml_contents(r)
    elif file_type == "fasta":
        _check_fasta_contents(r)

    with open(unp_file, "wb") as f:
        f.write(r.content)

def _check_json_response(response: requests.Response) -> None:
    """Checks if the response is a valid JSON.

    Args:
        response (requests.Response): Response to check.

    Raises:
        FileNotFoundError: If the response is not a valid JSON or if the entry is deleted.
    """
    content = response.json()
    try:
        if content["entryType"] == "Inactive":
            raise ObsoleteUniProtError(f"Uniprot entry {response.url} is deleted")
    except KeyError:
        raise AccessionNotFound(f"{response.url} response is invalid") from None


def _check_xml_contents(response: requests.Response) -> None:
    """Checks if the response is a valid XML by size."""

    root = ElementTree.fromstring(response.text)
    children = root.find(".//{http://uniprot.org/uniprot}entry")
    if children is None:
        raise ObsoleteUniProtError(
            f"Uniprot XML at {response.url} is probably blank."
            f"Size ({len(response.text)}) is too small."
            "Entry is probably deleted."
        )


def _check_fasta_contents(response: requests.Response) -> None:
    """Checks if the response is a valid FASTA by size."""
    if len(response.text) == 0:
        raise ObsoleteUniProtError(
            f"{response.url} is blank. Entry is probably deleted."
        )

def get_mismatches(seq1, seq2):
    return sum(aa1 != aa2 for aa1, aa2 in zip(seq1, seq2))

def get_identity(seq1, seq2):
    return round(sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2)) / float(len(seq2)), 2)

def get_coverage(query_from, query_to, query_len):
    return (query_to - query_from + 1) / query_len

def get_date():
    now = datetime.now()
    timestamp = now.strftime("%H_%d_%m_%Y")
    return timestamp

def make_path(base_dir: Path, id: str, sub_dir: str, filename: str, now: str = None) -> Path:
    """
    Creates a path with a timestamp, creates the parent folders,
    deletes the file if it exists, and returns the Path.
    """
    if now is None:
        now = get_date()

    dir_path = base_dir / f'{sub_dir}_{id}_{now}'
    full_path = dir_path / filename

    if dir_path.exists():
        shutil.rmtree(dir_path)

    dir_path.mkdir(parents=True, exist_ok=True)

    return full_path

def parse_extra_args(args: Optional[List[str]]) -> dict:
    kwargs = {}
    if not args:
        return kwargs

    i = 0
    while i < len(args):
        if args[i].startswith("-"):
            key = args[i].lstrip("-").replace("-", "_")
            # Handle flags without value (i.e., boolean flags)
            if i + 1 >= len(args) or args[i + 1].startswith("-"):
                kwargs[key] = True
                i += 1
            else:
                kwargs[key] = args[i + 1]
                i += 2
        else:
            i += 1
    return kwargs
