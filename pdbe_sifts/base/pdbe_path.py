"""Utility methods for getting paths to various files and directories to avoid hardcoding paths."""

from pathlib import Path
from pdbe_sifts.config import Config

conf = Config()


class SubDir:
    """Subdirectories in entry root for various itemsself.

    This should prevent hardcoding these in scripts.
    Can be accessed using dot notation and is mostly meant for use together with
    EntryLocation().get_entry_dir()
    """

    ARCHIVE_CIF_FTP = "mmCIF"
    CLEAN_CIF = "clean_mmcif"
    BCIF = "bcif"
    ARCHIVE_COPIES = "archive_copies"


def get_entry_dir(entry_id, base_dir, suffix=""):
    """The base directory for an entry where subfolders for specific items should exist."""
    return str(Path(base_dir, entry_id[1:3], entry_id, suffix))


def get_uniprot_cache_dir(unp_accession, base_dir=conf.cache.uniprot):
    """Path to UNP cache directory for specified UNP accession."""
    return str(Path(base_dir, unp_accession[0], unp_accession))


def get_entry_split_dir(entry_id, base_dir):
    """The base directory for an entry where subfolders for specific items should exist.

    This is used for things that take the same format as the archive divided directory.
    e.g. For 1atp, this would be $base_dir/at. Files can then be stored in here
    e.g. $base_dir/at/1atp.cif.gz
    """
    return str(Path(base_dir, entry_id[1:3]))


def get_clean_mmcif(entry_id, base_dir):
    """Path to Gzipped clean mmCIF file."""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.CLEAN_CIF,
            f"{entry_id}_new.cif.gz",
        )
    )

def get_updated_cif(entry_id, base_dir):
    """Path to updated mmCIF file produced at PDBe for specified entry."""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.CLEAN_CIF,
            f"{entry_id}_updated.cif.gz",
        )
    )


def get_archive_mmcif(entry_id, base_dir):
    """Path to work area copy of archive mmCIF file."""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir, suffix=SubDir.ARCHIVE_COPIES),
            f"{entry_id}.cif.gz",
        )
    )

def get_bcif(entry_id, base_dir):
    """Path to Binary mmCIF file for specified entry.

    Arguments:
        entry_id {str} -- PDB entry code
        base_dir {str} -- Base directory

    Returns:
        str -- Path to .bcif file. Does not check if file exists
    """
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.BCIF,
            f"{entry_id}.bcif",
        )
    )


def get_ftp_mmcif(entry_id, base_dir):
    """Get path to mmCIF-format entry file in FTP mirror."""
    return str(
        Path(
            base_dir,
            "data/structures/divided",
            SubDir.ARCHIVE_CIF_FTP,
            entry_id[1:3],
            f"{entry_id}.cif.gz",
        )
    )


def get_chem_comp_cif(het_code, base_dir, cvs=False):
    """Path to chemical component CIF file for specified hetcode."""
    if len(het_code) > 3 and cvs:
        return str(Path(base_dir, het_code[-2:], het_code, f"{het_code}.cif"))
    return str(Path(base_dir, het_code[0], het_code, f"{het_code}.cif"))


def get_sifts_res_csv(entry_id, base_dir):
    """Path to isoform residues CSV"""
    return str(
        Path(get_entry_dir(entry_id, base_dir, "sifts"), f"{entry_id}_res.csv.gz")
    )


def get_sifts_seg_csv(entry_id, base_dir):
    """Path to isoforms segments CSV"""
    return str(
        Path(get_entry_dir(entry_id, base_dir, "sifts"), f"{entry_id}_seg.csv.gz")
    )


def get_sifts_nf90_csv(entry_id, base_dir):
    """Path to NF90 segment CSV"""
    return str(
        Path(get_entry_dir(entry_id, base_dir, "sifts"), f"{entry_id}_nf90_seg.csv.gz")
    )


def get_sifts_updated_cif(entry_id, base_dir):
    """Path to mmCIF updated with SIFTS data"""
    return str(
        Path(get_entry_dir(entry_id, base_dir, "sifts"), f"{entry_id}_updated.cif.gz")
    )