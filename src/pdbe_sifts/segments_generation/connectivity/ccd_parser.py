#!/usr/bin/env python3

import os

import gemmi

from pdbe_sifts.base.paths import ccd_cache_path
from pdbe_sifts.base.utils import download_file_from_url

URL = "https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/ccd/"


def get_ccd_file(threeL_res: str) -> str:
    """Return the local path to the CCD CIF file for a residue, downloading it if needed.

    Args:
        threeL_res: Three-letter residue code (e.g. ``"ATP"``).

    Returns:
        Absolute path to the cached CCD ``.cif`` file.
    """
    threeL_res = threeL_res.upper()

    local_path = ccd_cache_path(threeL_res)

    if os.path.exists(local_path):
        return local_path

    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    url = URL + threeL_res[0] + "/" + threeL_res + "/" + threeL_res + ".cif"

    tmp_path = local_path + ".tmp"

    # Only download if another process isn't already doing it
    try:
        download_file_from_url(url, tmp_path)
        os.replace(tmp_path, local_path)  # atomic move
    except FileExistsError:
        # another process already finished the download
        pass
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    return local_path


def ccd_dict_builder(ccd_ters: dict) -> dict | int:
    """Build a compact terminal-atom descriptor from extracted CCD terminal data.

    Extracts the atom identifier and ordinal for both the N-terminal and
    C-terminal atoms from the dict returned by :meth:`CcdFile.extract_ters`.

    Args:
        ccd_ters: Dict with keys ``"n_ter"`` and/or ``"c_ter"``, each
            containing a ``(atom_index, atom_description_dict)`` tuple as
            returned by :meth:`CcdFile.extract_ters`.

    Returns:
        A dict ``{"N": (atom_id, pdbx_ordinal), "C": (atom_id, pdbx_ordinal)}``
        when both terminal atoms are present, or ``0`` if a required key is
        missing.
    """
    try:
        return {
            "N": (
                ccd_ters["n_ter"][1]["atom_id"],
                ccd_ters["n_ter"][1]["pdbx_ordinal"],
            ),
            "C": (
                ccd_ters["c_ter"][1]["atom_id"],
                ccd_ters["c_ter"][1]["pdbx_ordinal"],
            ),
        }
    except KeyError:
        return 0


class CcdFile:
    """Parse a CCD (Chemical Component Dictionary) CIF file for one residue.

    Provides access to N-terminal and C-terminal atom identifiers needed
    for the connectivity check step.
    """

    def __init__(self, threeL_res: str) -> None:
        """Initialise by fetching (or loading from cache) the CCD file.

        Args:
            threeL_res: Three-letter residue code (e.g. ``"SER"``).
        """
        self.threeL_res = threeL_res
        self.input_file = get_ccd_file(threeL_res)
        self.atoms = None

    def parse_ccd(self) -> None:
        """Parse the CCD CIF file and populate ``self.atoms`` with atom data."""
        cif = gemmi.cif.read_file(self.input_file)
        cif_block = cif.sole_block()
        cif_loop = cif_block.find_loop_item("_chem_comp_atom.comp_id")
        to_extract = [
            "_chem_comp_atom.comp_id",
            "_chem_comp_atom.atom_id",
            "_chem_comp_atom.alt_atom_id",
            "_chem_comp_atom.type_symbol",
            "_chem_comp_atom.pdbx_backbone_atom_flag",
            "_chem_comp_atom.pdbx_n_terminal_atom_flag",
            "_chem_comp_atom.pdbx_c_terminal_atom_flag",
            "_chem_comp_atom.pdbx_component_atom_id",
            "_chem_comp_atom.pdbx_component_comp_id",
            "_chem_comp_atom.pdbx_ordinal",
        ]

        columns = {col: cif_loop.loop.tags.index(col) for col in to_extract}
        atoms = {}
        for i in range(cif_loop.loop.length()):
            atom = {
                "comp_id": cif_loop.loop[i, columns["_chem_comp_atom.comp_id"]],
                "atom_id": cif_loop.loop[i, columns["_chem_comp_atom.atom_id"]],
                "alt_atom_id": cif_loop.loop[
                    i, columns["_chem_comp_atom.alt_atom_id"]
                ],
                "type_symbol": cif_loop.loop[
                    i, columns["_chem_comp_atom.type_symbol"]
                ],
                "pdbx_backbone_atom_flag": cif_loop.loop[
                    i, columns["_chem_comp_atom.pdbx_backbone_atom_flag"]
                ],
                "pdbx_n_terminal_atom_flag": cif_loop.loop[
                    i, columns["_chem_comp_atom.pdbx_n_terminal_atom_flag"]
                ],
                "pdbx_c_terminal_atom_flag": cif_loop.loop[
                    i, columns["_chem_comp_atom.pdbx_c_terminal_atom_flag"]
                ],
                "pdbx_component_atom_id": cif_loop.loop[
                    i, columns["_chem_comp_atom.pdbx_component_atom_id"]
                ],
                "pdbx_component_comp_id": cif_loop.loop[
                    i, columns["_chem_comp_atom.pdbx_component_comp_id"]
                ],
                "pdbx_ordinal": cif_loop.loop[
                    i, columns["_chem_comp_atom.pdbx_ordinal"]
                ],
            }
            atoms[i] = atom
        self.atoms = atoms

    def extract_ters(self) -> dict:
        """Extract the N-terminal and C-terminal atom descriptors from parsed atoms.

        Returns:
            Dict with keys ``"n_ter"`` and/or ``"c_ter"``, each mapping to
            ``(atom_index, atom_description_dict)``.
        """
        atoms_extracted = {}
        for atom, descrip in self.atoms.items():
            if (
                descrip["pdbx_n_terminal_atom_flag"] == "Y"
                and descrip["type_symbol"] == "N"
            ):
                atoms_extracted["n_ter"] = (atom, descrip)
            if (
                descrip["pdbx_c_terminal_atom_flag"] == "Y"
                and descrip["type_symbol"] == "C"
            ):
                atoms_extracted["c_ter"] = (atom, descrip)
        return atoms_extracted

    def process(self) -> dict | int:
        """Parse the CCD file and return terminal atom info.

        Returns:
            Dict mapping ``"N"`` and ``"C"`` to ``(atom_id, ordinal)`` tuples,
            or ``0`` if the required terminal atom fields are missing.
        """
        self.parse_ccd()
        return ccd_dict_builder(self.extract_ters())
