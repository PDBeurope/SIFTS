"""Align the canonical deposited sequence against the coordinate sequence.

Reads two sequences from an mmCIF file:

* **canonical** — the full sequence as deposited by the authors, from
  ``_entity_poly_seq``.
* **coordinate** — the sequence actually observed in the atomic coordinates,
  reconstructed from ``_atom_site`` ATOM records (model 1 only, unique
  residues ordered by auth_seq_id).

The two sequences are aligned with ``lalign36`` (local pairwise alignment)
and the result is returned as a structured dict containing both raw sequences,
the best alignment object, a human-readable annotation block, sequence
identity, and coordinate coverage.

Typical use-case: diagnosing missing loops, expression-tag residues, or
disordered regions that are absent from the electron-density map.
"""

import argparse
import logging
from pathlib import Path

from gemmi import cif

from pdbe_sifts.mmcif import _build_sequence
from pdbe_sifts.mmcif.chem_comp import ChemCompMapping
from pdbe_sifts.segments_generation.alignment import (
    annotate_alignment,
    do_alignment_lalign36,
    get_identity,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Low-level readers
# ---------------------------------------------------------------------------


def get_canonical_sequence(block, entity_id: str, cc: ChemCompMapping) -> str:
    """Return the canonical sequence for *entity_id* from ``_entity_poly_seq``.

    Args:
        block: A :class:`gemmi.cif.Block` object.
        entity_id: Entity identifier string (e.g. ``"1"``).
        cc: Chemical-component mapping used to convert three-letter codes to
            one-letter codes.

    Returns:
        One-letter amino acid sequence string.  Returns an empty string when
        the category is absent or the entity is not found.
    """
    entity_poly_seq = block.get_mmcif_category("_entity_poly_seq")
    if not entity_poly_seq:
        return ""
    return _build_sequence(entity_poly_seq, entity_id, cc)


def get_coordinate_sequence(
    block, entity_id: str, chain_id: str, cc: ChemCompMapping
) -> str:
    """Return the observed sequence for *entity_id* / *chain_id* from ``_atom_site``.

    Reads ATOM records from model 1, collects unique residues ordered by
    ``auth_seq_id`` + ``pdbx_PDB_ins_code``, and converts each three-letter
    ``label_comp_id`` to a one-letter code.  HETATM records are excluded.

    Args:
        block: A :class:`gemmi.cif.Block` object.
        entity_id: Entity identifier string (e.g. ``"1"``).
        chain_id: Author asymmetric-unit chain identifier (e.g. ``"A"``).
        cc: Chemical-component mapping used to convert three-letter codes to
            one-letter codes.

    Returns:
        One-letter amino acid sequence string of the observed residues.
        Returns an empty string when no matching ATOM records are found.
    """
    atom_site = block.get_mmcif_category("_atom_site")
    if not atom_site:
        return ""

    groups = atom_site.get("group_PDB", [])
    entity_ids = atom_site.get("label_entity_id", [])
    chain_ids = atom_site.get("auth_asym_id", [])
    seq_ids = atom_site.get("auth_seq_id", [])
    ins_codes = atom_site.get("pdbx_PDB_ins_code", [])
    comp_ids = atom_site.get("label_comp_id", [])
    model_nums = atom_site.get("pdbx_PDB_model_num", [])

    seen: dict[tuple, str] = {}
    for i, group in enumerate(groups):
        if group != "ATOM":
            continue
        if entity_ids[i] != entity_id:
            continue
        if chain_ids[i] != chain_id:
            continue
        # model 1 only
        if model_nums and model_nums[i] not in ("1", 1):
            continue

        ins = ins_codes[i] if ins_codes else "."
        key = (int(seq_ids[i]), ins)
        if key not in seen:
            seen[key] = comp_ids[i]

    if not seen:
        return ""

    return "".join(cc.get(seen[k]) for k in sorted(seen))


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------


class Seq2Seq:
    """Align the canonical deposited sequence against the coordinate sequence.

    Reads both sequences from a single mmCIF file and runs a local pairwise
    alignment (``lalign36``) between them.

    Args:
        cif_file: Path to the mmCIF file (plain or gzip-compressed).
        entity_id: Entity identifier, as a string (e.g. ``"1"``).
        chain_id: Author asymmetric-unit chain identifier (e.g. ``"A"``).

    Attributes:
        cif_file: Resolved :class:`pathlib.Path` to the input CIF.
        entity_id: Entity identifier string.
        chain_id: Chain identifier string.
    """

    def __init__(self, cif_file, entity_id: str, chain_id: str) -> None:
        """Initialise Seq2Seq with a CIF file and target entity/chain.

        Args:
            cif_file: Path to the mmCIF file.
            entity_id: Entity identifier string (e.g. ``"1"``).
            chain_id: Author chain identifier (e.g. ``"A"``).
        """
        self.cif_file = Path(cif_file)
        self.entity_id = str(entity_id)
        self.chain_id = str(chain_id)

    def run(self) -> dict:
        """Execute the canonical-vs-coordinate alignment.

        Reads both sequences from the CIF file, aligns them with
        ``lalign36``, and returns a summary dict.

        Returns:
            A dict with the following keys:

            * ``canonical``  (:class:`str`) — full deposited sequence
              from ``_entity_poly_seq``.
            * ``coordinate`` (:class:`str`) — observed sequence from
              ``_atom_site`` ATOM records.
            * ``alignment``  (:class:`Bio.Align.MultipleSeqAlignment` or
              ``None``) — best (first) alignment object returned by
              lalign36; ``None`` if alignment failed.
            * ``annotated``  (:class:`str`) — human-readable three-line
              UNP/ALG/PDB block; empty string on failure.
            * ``identity``   (:class:`float`) — fraction of identical
              positions in the best alignment (0.0–1.0).
            * ``coverage``   (:class:`float`) — fraction of the canonical
              sequence covered by the aligned coordinate sequence (0.0–1.0).

        Raises:
            FileNotFoundError: If :attr:`cif_file` does not exist.
        """
        if not self.cif_file.exists():
            raise FileNotFoundError(f"CIF file not found: {self.cif_file}")

        cc = ChemCompMapping()
        block = cif.read(str(self.cif_file)).sole_block()

        canonical = get_canonical_sequence(block, self.entity_id, cc)
        coordinate = get_coordinate_sequence(
            block, self.entity_id, self.chain_id, cc
        )

        if not canonical:
            logger.warning(
                f"No canonical sequence found for entity {self.entity_id}"
            )
        if not coordinate:
            logger.warning(
                f"No coordinate sequence found for entity {self.entity_id} "
                f"chain {self.chain_id}"
            )

        if not canonical or not coordinate:
            return {
                "canonical": canonical,
                "coordinate": coordinate,
                "alignment": None,
                "annotated": "",
                "identity": 0.0,
                "coverage": 0.0,
            }

        try:
            alignments = list(do_alignment_lalign36(canonical, coordinate))
            best = alignments[0] if alignments else None
        except Exception as exc:
            logger.error(f"lalign36 failed: {exc}")
            best = None

        if best is None:
            return {
                "canonical": canonical,
                "coordinate": coordinate,
                "alignment": best,
                "annotated": "",
                "identity": 0.0,
                "coverage": 0.0,
            }

        al_canonical = str(best[0].seq)  # aligned canonical  (with gaps)
        al_coordinate = str(best[1].seq)  # aligned coordinate (with gaps)

        identity = get_identity(al_canonical, al_coordinate)

        # coverage = aligned coordinate length (no gaps) / canonical length
        coord_aligned_len = len(al_coordinate.replace("-", ""))
        coverage = (
            round(coord_aligned_len / len(canonical), 2) if canonical else 0.0
        )

        annotated = annotate_alignment(al_canonical, al_coordinate)

        return {
            "canonical": canonical,
            "coordinate": coordinate,
            "alignment": best,
            "annotated": annotated,
            "identity": identity,
            "coverage": coverage,
        }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def run() -> None:
    """Command-line entry point: ``pdbe_sifts seq2seq``.

    Args (parsed from command line):
        -i / --input-cif:  Path to the mmCIF file.
        -e / --entity-id:  Entity identifier (e.g. ``1``).
        -c / --chain-id:   Author chain identifier (e.g. ``A``).
    """
    parser = argparse.ArgumentParser(
        prog="pdbe_sifts seq2seq",
        description="Align the canonical deposited sequence against the coordinate sequence.",
    )
    parser.add_argument(
        "-i", "--input-cif", required=True, help="mmCIF file path"
    )
    parser.add_argument(
        "-e", "--entity-id", required=True, help="Entity ID (e.g. 1)"
    )
    parser.add_argument(
        "-c", "--chain-id", required=True, help="Author chain ID (e.g. A)"
    )
    args = parser.parse_args()

    result = Seq2Seq(args.input_cif, args.entity_id, args.chain_id).run()

    print(f"Entity {args.entity_id}  chain {args.chain_id}")
    print(f"Canonical length  : {len(result['canonical'])}")
    print(f"Coordinate length : {len(result['coordinate'])}")
    print(f"Identity          : {result['identity']:.2f}")
    print(f"Coverage          : {result['coverage']:.2f}")
    print()
    if result["annotated"]:
        print(result["annotated"])
    else:
        print("No alignment produced.")
