"""
Author: Preeti Choudhary
Desription: Parse the ccd cif file and generate three_to_one_letter_mapping.csv file
"""

import csv
import tempfile
from datetime import datetime
from email.utils import parsedate_to_datetime
from pathlib import Path

import gemmi
import requests

from pdbe_sifts.base import utils


def parse_ccd_file(ccd_cif_file, ccd_mapping_file):
    """Parse the PDBe CCD cif file and generate three to one letter mapping file.

    Args:
        ccd_cif_file: Path to CCD CIF file (components.cif from PDBe FTP).
        ccd_mapping_file: Output path for three_to_one_letter_mapping.csv.
    """
    # Open output three_to_one_letter_mapping.csv file for writing
    print(f"Writing the three_to_one_mapping.csv: {ccd_mapping_file}")
    with open(ccd_mapping_file, "w", newline="") as f:
        writer = csv.writer(f)
        # Read PDBe CCD CIF file
        # make sure file path is a string as gemmi does not accept a Path object
        doc = gemmi.cif.read_file(str(ccd_cif_file))

        for block in doc:
            # get the chemp_comp category
            chem_comp_cat = block.get_mmcif_category("_chem_comp")

            if chem_comp_cat:
                # get the one/three letter code for a given ccd
                one_letter_list = chem_comp_cat.get("one_letter_code")
                three_letter_list = chem_comp_cat.get("three_letter_code")

                if one_letter_list and three_letter_list:
                    one = one_letter_list[0]
                    three = three_letter_list[0]

                    # If missing values i.e ("?" or ".") use X as one letter code
                    if one is None:
                        one = "X"
                    writer.writerow([three, one])


# Path of the CCD cif file from PDBe FTP area
DEFAULT_CCD_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/ccd/components.cif"
)


def get_remote_ccd_date(url: str = DEFAULT_CCD_URL) -> datetime:
    """Return the Last-Modified datetime of the remote CCD file (UTC-aware).

    Raises:
        ValueError: If the server does not return a Last-Modified header.
        requests.HTTPError: If the HEAD request fails.
    """
    response = requests.head(url, allow_redirects=True, timeout=30)
    response.raise_for_status()
    last_modified = response.headers.get("Last-Modified")
    if not last_modified:
        raise ValueError(f"No Last-Modified header returned from {url}")
    return parsedate_to_datetime(last_modified)


def get_cached_ccd_date(readme_path: Path) -> datetime | None:
    """Read the generation date from the README sidecar file.

    Returns None if the file doesn't exist or the date line is missing/malformed.
    """
    if not readme_path.exists():
        return None
    try:
        for line in readme_path.read_text().splitlines():
            if line.startswith("Generated:"):
                date_str = line.split(":", 1)[1].strip()
                return datetime.fromisoformat(date_str)
    except Exception:
        return None
    return None


def generate_mapping_to_cache(
    url: str = DEFAULT_CCD_URL, csv_path: Path | None = None
):
    """Download the CCD CIF file and regenerate the three-to-one mapping CSV.

    Writes the CSV to *csv_path* and a README sidecar recording the source date.
    The README path is derived from the CSV path by replacing the extension with
    ``.README``.

    Args:
        url: URL of the CCD CIF file.
        csv_path: Destination path for the CSV. When None, the configured cache
            path is used (``conf.cache.three_to_one``).

    Raises:
        ValueError: If no csv_path is provided and the config key is unset.
    """
    if csv_path is None:
        from pdbe_sifts.base.paths import get_conf_three_to_one_csv_path

        csv_path = get_conf_three_to_one_csv_path()
        if csv_path is None:
            raise ValueError(
                "conf.cache.three_to_one is not set and no csv_path was provided."
            )

    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    readme_path = csv_path.with_suffix(".README")

    with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        utils.download_file_from_url(url, tmp_path)
        parse_ccd_file(tmp_path, csv_path)
    finally:
        tmp_path.unlink(missing_ok=True)

    remote_date = get_remote_ccd_date(url)
    readme_path.write_text(
        f"Generated: {remote_date.isoformat()}\nSource: {url}\n"
    )


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate three_to_one_letter_mapping.csv from CCD CIF file"
    )

    parser.add_argument(
        "--ccd_file",
        help="Path to CCD CIF file (if not provided, it will be downloaded)",
    )

    parser.add_argument(
        "--output_path",
        required=True,
        help="Path to output three_to_one_mapping.csv file",
    )

    parser.add_argument(
        "--ccd_url",
        default=DEFAULT_CCD_URL,
        help="URL to download CCD file if --ccd_file is not provided",
    )

    args = parser.parse_args()
    # Decide where CCD file comes from
    if args.ccd_file:
        ccd_file = Path(args.ccd_file)
    else:
        # Default download location - current directory from where code is being run
        base_dir = Path.cwd()
        ccd_file = Path(base_dir, "components.cif")
        if not ccd_file.exists():
            utils.download_file_from_url(args.ccd_url, ccd_file)
        else:
            print(f"Using existing ccd file: {ccd_file}")

    output_path = Path(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = Path(output_path, "three_to_one_letter_mapping.csv")

    parse_ccd_file(ccd_file, output_file)


if __name__ == "__main__":
    main()
