#!/usr/bin/env python3

import argparse
import os
from datetime import date, datetime
from pathlib import Path
from xml.etree import ElementTree

import funcy
import requests
from dateutil.relativedelta import WE, relativedelta

from pdbe_sifts.base.exceptions import AccessionNotFound, ObsoleteUniProtError
from pdbe_sifts.base.log import logger
from pdbe_sifts.base.paths import uniprot_cache_dir as get_uniprot_cache_dir

UNIPROT_REGEX = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
UNIPROT_API_BASE_URL = "https://rest.uniprot.org/uniprotkb"


def fetch_uniprot_file(uniprot_id: str, file_type: str, fail_silently=False) -> str:
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
        unp_dir = get_uniprot_cache_dir(uniprot_id)
        if file_type not in ["xml", "json", "fasta"]:
            raise ValueError(f"Invalid file type {file_type}. Must be one of xml, json, fasta")
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

    Raises:
        ObsoleteUniProtError: If the entry is deleted.
        AccessionNotFound: If the response is invalid.
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
        raise ObsoleteUniProtError(f"{response.url} is blank. Entry is probably deleted.")


def get_date():
    now = datetime.now()
    timestamp = now.strftime("%H_%d_%m_%Y")
    return timestamp


def get_next_release_date() -> date:
    """Returns the next PDBe release date.

    This should be the next Wednesday after the today's date unless today is Wednesday
    in which case it returns same as `date.today()`.
    """
    return date.today() + relativedelta(weekday=WE(+1))


def get_allocated_cpus():
    """Return the number of CPUs available to this process.

    Environment variables (checked in order):
        SLURM_CPUS_PER_TASK — set automatically by SLURM; used on HPC clusters.
        Falls back to os.sched_getaffinity (Linux) or os.cpu_count().
    """
    if "SLURM_CPUS_PER_TASK" in os.environ:
        return int(os.environ["SLURM_CPUS_PER_TASK"])
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def get_cpu_count():
    """Return the number of parallel threads to use for internal alignment jobs.

    Environment variables (checked in order):
        SIFTS_N_PROC — override the CPU count (set by segments_batch per worker).
        Falls back to get_allocated_cpus().
    """
    if "SIFTS_N_PROC" in os.environ:
        return int(os.environ["SIFTS_N_PROC"])
    return get_allocated_cpus()


class SiftsAction(argparse.Action):
    def __init__(self, envvar=None, confvar=None, required=False, default=None, **kwargs):
        if not default:
            if envvar and envvar in os.environ:
                default = os.environ[envvar]
            if not default and confvar:
                default = confvar

        if required and default:
            required = False
        super().__init__(default=default, required=required, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)


def download_file_from_url(url: str, output_file: str) -> None:
    """Downloads a file from a URL.

    Args:
        url (str): URL to download from.
        output_file (str): Path to output file.
    """
    logger.info(f"Downloading {url} to {output_file}")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(output_file, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
