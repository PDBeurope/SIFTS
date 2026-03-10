"""Utility methods for getting paths to various files and directories to avoid hardcoding paths."""

import re
from pathlib import Path
from xml.etree import ElementTree as ET

from pdbe_sifts.config import load_config

conf = load_config()


def is_valid_emdb_code(emdb_id):
    return re.match(r"EMD[-_]\d{4,5}", emdb_id, re.IGNORECASE)


class SubDir:
    """Subdirectories in entry root for various itemsself.

    This should prevent hardcoding these in scripts.
    Can be accessed using dot notation and is mostly meant for use together with
    EntryLocation().get_entry_dir()

    Examples:
        >>> from segments_gen.base.pdbe_path import SubDir, EntryLocation
        >>> import os
        >>> SubDir.IMAGES
        'images'
        >>> # Get images directory for an entry
        >>> os.path.join(EntryLocation().get_entry_dir('1atp'), SubDir.IMAGES )
        '/nfs/msd/release/data/entry/at/1atp/images'

    """

    IMAGES = "images"
    ASSEMBLIES = "assembly_generation"
    ARCHIVE_ASSEMBLIES = "assemblies"
    STRUCTURE_FACTORS = "structure_factors"
    PDB = "pdb"
    ARCHIVE_CIF_FTP = "mmCIF"
    ARCHIVE_CIF_COPY = "original_mmcif"
    CLEAN_CIF = "clean_mmcif"
    INTERACTIONS = "interactions"
    CITATION_MAPPING = "citation_mapping_and_images"
    BINDING_SITE = "binding_site_update"
    DOSS = "doss"
    EDS = "eds"
    TOPOLOGY = "topology"
    ASSEMBLY_SYMMETRY = "assembly_symmetry"
    PISA = "pisa"
    BCIF = "bcif"
    CITATION_UPDATE = "citation_update"
    ARCHIVE_COPIES = "archive_copies"
    DSSP = "dssp"
    BOUNDMOLECULES = "bound_molecules"
    MOLE = "mole"


archive_divided = "data/structures/divided"


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


def get_dssp_cif_file(entry_id, base_dir):
    """Path to dssp mmCIF file."""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.DSSP,
            f"{entry_id}.cif",
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


def get_archive_pdb(entry_id, base_dir):
    """Path to work area copy of wwPDB ent file."""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir, suffix=SubDir.ARCHIVE_COPIES),
            f"pdb{entry_id}.ent.gz",
        )
    )


def get_archive_structure_factor_file(entry_id, base_dir):
    """Get path to structure factor file in archive."""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir, suffix=SubDir.ARCHIVE_COPIES),
            f"r{entry_id}sf.ent.gz",
        )
    )


def get_archive_validation_xml(entry_id, base_dir):
    """Path to work area copy of wwPDB validation XML file."""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir, suffix=SubDir.ARCHIVE_COPIES),
            f"{entry_id}_validation.xml.gz",
        )
    )


def get_archive_assembly_file(entry_id, assembly_id, base_dir=None):
    """Get path to working area copy of FTP assembly file.

    Args:
        entry_id(str): PDB entry code
        assembly_id(int): Assembly ID
        base_dir: Entry base directory

    Returns:
        pathlib.Path: Path to assembly file
    """
    return Path(
        get_entry_dir(entry_id, base_dir, suffix=SubDir.ARCHIVE_ASSEMBLIES),
        f"{entry_id}-assembly{assembly_id}.cif.gz",
    )


def get_archive_assembly_files(entry_id, base_dir=None) -> list[tuple[int, str]]:
    """Get paths to all copies of FTP assembly files in work area.

    Args:
        entry_id(str): PDB entry code
        base_dir: Entry base directory

    Returns:
        pathlib.Path: Path to assembly file
    """
    files = Path(
        get_entry_dir(entry_id, base_dir, suffix=SubDir.ARCHIVE_ASSEMBLIES)
    ).glob(f"{entry_id}-assembly*.cif.gz")
    output = []
    for filepath in files:
        match = re.match(rf"{entry_id}-assembly(\d+).cif.gz", filepath.name)
        if match:
            assembly_id = int(match.group(1))
            output.append((assembly_id, str(filepath)))

    return output


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
            archive_divided,
            SubDir.ARCHIVE_CIF_FTP,
            entry_id[1:3],
            f"{entry_id}.cif.gz",
        )
    )


def get_ftp_pdb_ent(entry_id, base_dir):
    """Get path to PDB-format entry file in FTP mirror.

    Args:
        entry_id(str): PDB entry code
        base_dir: Root of ftp mirror typically a pdb folder
    """
    return str(
        Path(
            base_dir,
            archive_divided,
            SubDir.PDB,
            entry_id[1:3],
            f"pdb{entry_id}.ent.gz",
        )
    )


def get_ftp_validation_xml(entry_id, base_dir):
    """Gets path to validation XML file in FTP mirror.

    Currently in $FTPDIR/pdb/validation_reports/at/1atp/1atp_validation.xml.gz
    Args:
        entry_id(str): PDB entry code
        base_dir: Root of ftp mirror typically a pdb folder
    """
    return str(
        Path(
            get_entry_dir(entry_id, Path(base_dir, "validation_reports")),
            f"{entry_id}_validation.xml.gz",
        )
    )


def get_ftp_structure_factor_file(entry_id, base_dir):
    """Get path to structure factor file in FTP mirror.

    Args:
        entry_id(str): PDB entry code
        base_dir: Root of ftp mirror typically a pdb folder
    """
    return str(
        Path(
            base_dir,
            archive_divided,
            SubDir.STRUCTURE_FACTORS,
            entry_id[1:3],
            f"r{entry_id}sf.ent.gz",
        )
    )


def get_chem_comp_cif(het_code, base_dir, cvs=False):
    """Path to chemical component CIF file for specified hetcode."""
    if len(het_code) > 3 and cvs:
        return str(Path(base_dir, het_code[-2:], het_code, f"{het_code}.cif"))
    return str(Path(base_dir, het_code[0], het_code, f"{het_code}.cif"))


def get_superpose_dir(unp_accession, base_dir):
    """Path to superpose output for uniprot."""
    return str(Path(base_dir, unp_accession[0], unp_accession))


def get_prd_cif(prdcc_code, base_dir):
    """Path to PRD CIF file from wwPDB for specified PRD/PRDCC code."""
    prefix, code = prdcc_code.split("_")
    if prefix == "PRD":
        prdcc_code = f"{prefix}CC_{code}"
    one_letter_code = prdcc_code[-1]
    return str(Path(base_dir, one_letter_code, f"{prdcc_code}.cif"))


def get_prd_cif_enriched(prdcc_code, base_dir):
    """Path to enriched PRD CIF file for specified PRD/PRDCC code."""
    prefix, code = prdcc_code.split("_")
    if prefix == "PRD":
        prdcc_code = f"{prefix}CC_{code}"
    one_letter_code = prdcc_code[-1]
    return str(Path(base_dir, one_letter_code, prdcc_code, f"{prdcc_code}.cif"))


def get_clc_dir(clc_id, base_dir):
    """ "Path to CLC directory for specified CLC ID"""
    return str(Path(base_dir, clc_id[-1], clc_id))


def get_clc_cif(clc_id, base_dir):
    """Path to CLC CIF file for specified CLC ID."""
    clc_dir = get_clc_dir(clc_id, base_dir)
    return str(Path(clc_dir, f"{clc_id}.cif"))


def get_clc_reader_result(clc_id, base_dir):
    """Return clc_reader_result for specified CLC ID"""
    clc_dir = get_clc_dir(clc_id, base_dir)
    return str(Path(clc_dir, "clc_reader_results.pkl"))


def get_assembly_dir(entry_id, base_dir=None):
    """Path to assembly directory for specified entry."""
    return str(Path(get_entry_dir(entry_id, base_dir), SubDir.ARCHIVE_ASSEMBLIES))


def get_assembly_files(entry_id, base_dir=None):
    """Get list of assembly file paths for specified entry."""
    assembly_dir = Path(get_assembly_dir(entry_id, base_dir))
    return list(assembly_dir.glob(f"{entry_id}-assembly*.cif.gz"))


def get_assembly_dssp_files(entry_id, base_dir=None):
    """Get list of assembly DSSP file paths for specified entry."""
    dssp_dir = Path(get_entry_dir(entry_id, base_dir), SubDir.DSSP)
    return list(dssp_dir.glob(f"{entry_id}_assembly*.cif"))


def get_assembly_xml(entry_id, base_dir=None):
    """Path to assembly XML for specified entry."""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.ARCHIVE_ASSEMBLIES,
            f"{entry_id}-assembly.xml",
        )
    )


def get_ftp_assembly_file(entry_id, assembly_id, base_dir=None):
    """Get path to assembly file in FTP mirror.

    Args:
        entry_id(str): PDB entry code
        assembly_id(int): Assembly ID
        base_dir: Root of ftp mirror typically a `pdb` folder

    Returns:
        pathlib.Path: Path to assembly file
    """
    return Path(
        base_dir,
        "data/assemblies/mmCIF/divided",
        entry_id[1:3],
        f"{entry_id}-assembly{assembly_id}.cif.gz",
    )


def get_ftp_assembly_files(entry_id, base_dir=None):
    """Get list of paths to assembly files in FTP mirror.

    Args:
        entry_id(str): PDB entry code
        base_dir: Root of ftp mirror typically a `pdb` folder

    Returns:
        list[pathlib.Path]: List of paths to assembly files
    """
    files = Path(base_dir, "data/assemblies/mmCIF/divided", entry_id[1:3]).glob(
        f"{entry_id}-assembly*.cif.gz"
    )
    output = []
    for filepath in files:
        match = re.match(rf"{entry_id}-assembly(\d+).cif.gz", filepath.name)
        if match:
            assembly_id = int(match.group(1))
            output.append((assembly_id, filepath))

    return output


def get_assemblies_info(entry_id, base_dir=None):
    """Get detailed assembly info for each of the generated assemblies.
    As of now these are the available attributes:
        * composition
        * id
        * molecular_weight
        * name
        * order
        * prefered (funny spelling)
        * type

    Args:
        entry_id (str): PDB entry id.
        base_dir (str, optional): Base directory to use instead of
            default $DATA_ENTRY_DIR. Defaults to None.

    Raises:
        AttributeError: If assembly path does not exist.

    Returns:
        list[dict[str,str]]: Assemblies meta information.
    """
    assembly_path = get_assembly_xml(entry_id, base_dir)

    if not Path(assembly_path).is_file():
        raise AttributeError(
            f"Assembly configuration file {assembly_path} does not exist."
        )

    xml = ET.parse(assembly_path)

    return [node.attrib for node in xml.getroot().iter("assembly")]


def get_preferred_assembly_info(entry_id, base_dir):
    """Shorthand to get assembly information in a dictionary form.

    Returns:
        dict[str,str]: Preferred assembly meta information.
                    keys: composition,id,moleculular_weight,name,order,prefered,type
    """
    assemblies = get_assemblies_info(entry_id, base_dir)

    return next((x for x in assemblies if x["prefered"] == "True"), {})


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


def get_eds_mtz(entry_id, base_dir):
    """Path to MTZ file from EDS"""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.EDS,
            f"{entry_id}_map.mtz",
        )
    )


def get_eds_xml(entry_id, base_dir):
    """Path to XML result file from EDS"""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.EDS,
            f"{entry_id}_eds.xml",
        )
    )


def get_eds_2fofc_ccp4_map(entry_id, base_dir):
    """Path to 2Fo-Fc CCP4 map file from EDS"""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.EDS,
            f"{entry_id}.ccp4",
        )
    )


def get_eds_fofc_ccp4_map(entry_id, base_dir):
    """Path to Fo-Fc CCP4 map file from EDS"""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.EDS,
            f"{entry_id}_diff.ccp4",
        )
    )


def get_volume_server_map(entry_id, base_dir):
    """Path to volume server map file"""
    return str(
        Path(
            get_entry_dir(entry_id, base_dir),
            SubDir.EDS,
            f"{entry_id}.mdb",
        )
    )


def get_pisa_xml_files(entry_id, assembly_id, base_dir):
    """Returns a tuple containing paths to assembly and interface XMLs for the assembly"""

    xml_dir = Path(get_entry_dir(entry_id, base_dir), SubDir.PISA, "xmls")
    return (
        str(xml_dir / f"{entry_id}_{assembly_id}_assembly.xml"),
        str(xml_dir / f"{entry_id}_{assembly_id}_interfaces.xml"),
    )


def get_pisa_json_files(entry_id, assembly_id, base_dir):
    """Returns a tuple containing paths to assembly and interface JSONs for the assembly"""

    pisa_dir = Path(get_entry_dir(entry_id, base_dir), SubDir.PISA)
    return (
        str(pisa_dir / f"{entry_id}_assembly{assembly_id}.json"),
        str(pisa_dir / f"{entry_id}_assembly{assembly_id}_interfaces.json"),
    )


def get_mole_json_file(entry_id, assembly_id, base_dir):
    """Returns a string containing the path to the Mole json for the specified assembly ID"""

    mole_dir = Path(get_entry_dir(entry_id, base_dir), SubDir.MOLE)
    return str(mole_dir / f"{entry_id}-assembly{assembly_id}.json")


def get_emdb_number(emdb_id):
    """Validates EMDB ID and returns the numerical part of the ID.

    Args:
        emdb_id (str): EMDB ID.

    Raises:
        ValueError: If EMDB ID is not valid.
    """
    if not is_valid_emdb_code(emdb_id):
        raise ValueError(
            f"Invalid EMDB ID: {emdb_id}. Format should be EMD-#### or EMD_####."
        )
    # Split on underscore or hyphen
    return re.split("_|-", emdb_id)[1]


def get_emdb_map(emdb_id, base_dir):
    """Path to EMDB map file"""
    emd_id = get_emdb_number(emdb_id)

    return str(Path(base_dir, f"EMD-{emd_id}", "map", f"emd_{emd_id}.map.gz"))


def get_em_volume_mdb(emdb_id, base_dir):
    """Path to Molstar density map file"""
    emd_id = get_emdb_number(emdb_id)
    return str(Path(base_dir, f"EMD-{emd_id}", f"emd-{emd_id}.mdb"))


def get_ftp_validation_2fofc_map_coeff(entry_id, base_dir):
    """Gets path to validation 2FO-FC map coeffient file in FTP mirror.

    Currently in $FTPDIR/pdb/validation_reports/at/1atp/1atp_validation_2fo-fc_map_coef.cif.gz
    Args:
        entry_id(str): PDB entry code
        base_dir: Root of ftp mirror typically a pdb folder
    """
    return str(
        Path(
            get_entry_dir(entry_id, Path(base_dir, "validation_reports")),
            f"{entry_id}_validation_2fo-fc_map_coef.cif.gz",
        )
    )


def get_ftp_validation_fofc_map_coeff(entry_id, base_dir):
    """Gets path to validation FO-FC map coeffient file in FTP mirror.

    Currently in $FTPDIR/pdb/validation_reports/at/1atp/1atp_validation_fo-fc_map_coef.cif.gz
    Args:
        entry_id(str): PDB entry code
        base_dir: Root of ftp mirror typically a pdb folder
    """
    return str(
        Path(
            get_entry_dir(entry_id, Path(base_dir, "validation_reports")),
            f"{entry_id}_validation_fo-fc_map_coef.cif.gz",
        )
    )
