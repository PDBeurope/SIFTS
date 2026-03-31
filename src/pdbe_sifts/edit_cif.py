"""Enrich an mmCIF file by injecting missing categories required by the pipeline.

Some CIF files — in-house models, partial depositions, or files exported by
third-party tools — lack categories that the SIFTS pipeline expects.  This
module provides a pre-processing step that reads the original CIF, computes
the missing data from what *is* present, and writes an enriched copy.

Currently supported additions:

* ``_pdbx_poly_seq_scheme`` — reconstructed from ``_entity_poly_seq`` and
  ``_atom_site`` via a local pairwise alignment (``lalign36``).

Usage::

    # Python API — process all polypeptide entity/chain pairs automatically
    from pdbe_sifts.edit_cif import EditCif

    ec = EditCif("1abc.cif", output_dir="./out/")
    ec.add_pdbx_poly_seq_scheme()
    path = ec.write()           # → ./out/1abc_edited.cif

    # Python API — process a single entity/chain pair
    ec = EditCif("1abc.cif", output_dir="./out/", entity_id="1", chain_id="A")
    ec.add_pdbx_poly_seq_scheme()
    path = ec.write()

    # CLI — all entities (default)
    pdbe_sifts edit_cif -i 1abc.cif -o ./out/

    # CLI — single entity/chain pair
    pdbe_sifts edit_cif -i 1abc.cif -e 1 -c A -o ./out/
"""

import logging
from pathlib import Path

from gemmi import cif

from pdbe_sifts.mmcif.chem_comp import ChemCompMapping

logger = logging.getLogger(__name__)


class EditCif:
    """Pre-process a CIF file by injecting missing mmCIF categories.

    Reads the CIF on construction, accumulates enrichments via ``add_*``
    methods, and writes the result with :meth:`write`.

    When *entity_id* and *chain_id* are both ``None`` (the default), all
    polypeptide entity/chain pairs found in the CIF are processed
    automatically.  Pass both to restrict reconstruction to a single pair.

    Args:
        cif_file: Path to the input mmCIF file (plain or gzip-compressed).
        output_dir: Directory where the enriched CIF will be written.
        entity_id: Target entity identifier string (e.g. ``"1"``).
            ``None`` means process all entities.
        chain_id: Author asymmetric-unit chain identifier (e.g. ``"A"``).
            ``None`` means process all chains.

    Attributes:
        cif_file: Resolved :class:`pathlib.Path` to the input CIF.
        entity_id: Entity identifier string, or ``None``.
        chain_id: Chain identifier string, or ``None``.
        output_dir: Resolved :class:`pathlib.Path` to the output directory.

    Raises:
        ValueError: If exactly one of *entity_id* / *chain_id* is provided
            (both or neither must be supplied).
    """

    def __init__(
        self,
        cif_file,
        output_dir,
        entity_id: str | None = None,
        chain_id: str | None = None,
    ) -> None:
        """Initialise EditCif and parse the input CIF block.

        Args:
            cif_file: Path to the input mmCIF file.
            output_dir: Output directory path.
            entity_id: Entity identifier string (e.g. ``"1"``), or ``None``
                to process all polypeptide entities.
            chain_id: Author chain identifier (e.g. ``"A"``), or ``None``
                to process all chains.

        Raises:
            FileNotFoundError: If *cif_file* does not exist.
            ValueError: If exactly one of *entity_id* / *chain_id* is given.
        """
        if (entity_id is None) != (chain_id is None):
            raise ValueError(
                "entity_id and chain_id must both be provided or both be omitted."
            )

        self.cif_file = Path(cif_file)
        self.entity_id = str(entity_id) if entity_id is not None else None
        self.chain_id = str(chain_id) if chain_id is not None else None
        self.output_dir = Path(output_dir)

        if not self.cif_file.exists():
            raise FileNotFoundError(f"CIF file not found: {self.cif_file}")

        self._doc = cif.read(str(self.cif_file))
        self._block = self._doc.sole_block()
        self._cc = ChemCompMapping()

    # ------------------------------------------------------------------
    # Enrichment methods
    # ------------------------------------------------------------------

    def add_pdbx_poly_seq_scheme(self) -> bool:
        """Reconstruct and inject ``_pdbx_poly_seq_scheme`` if absent.

        Uses ``_entity_poly_seq`` (with ``_entity_poly`` as fallback) and
        ``_atom_site`` as sources, aligned via ``lalign36``, to build the
        full residue-level scheme table.  Does nothing if the category is
        already present.

        When :attr:`entity_id` and :attr:`chain_id` are set, only that
        specific pair is processed.  Otherwise, all polypeptide entity/chain
        pairs found in the CIF are processed automatically.

        Returns:
            ``True`` when the category was added, ``False`` when it was
            already present or the reconstruction failed.
        """
        from pdbe_sifts.sifts_to_mmcif.poly_seq_scheme import (
            build_pdbx_poly_seq_scheme,
            build_pdbx_poly_seq_scheme_all,
            write_pdbx_poly_seq_scheme,
        )

        if self._block.get_mmcif_category("_pdbx_poly_seq_scheme"):
            logger.warning(
                "_pdbx_poly_seq_scheme already present — skipping reconstruction"
            )
            return False

        if self.entity_id and self.chain_id:
            logger.info(
                f"Reconstructing _pdbx_poly_seq_scheme for entity {self.entity_id} "
                f"chain {self.chain_id} …"
            )
            scheme = build_pdbx_poly_seq_scheme(
                self._block, self.entity_id, self.chain_id, self._cc
            )
        else:
            logger.info(
                "Reconstructing _pdbx_poly_seq_scheme for all polypeptide chains …"
            )
            scheme = build_pdbx_poly_seq_scheme_all(self._block, self._cc)

        if not scheme:
            logger.error("Reconstruction returned empty — nothing written")
            return False

        write_pdbx_poly_seq_scheme(self._block, scheme)
        logger.info(
            f"_pdbx_poly_seq_scheme injected ({len(scheme['seq_id'])} rows)"
        )
        return True

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def write(self) -> Path:
        """Write the enriched CIF to *output_dir*.

        The output filename is ``{entry_id}_edited.cif`` where *entry_id* is
        derived from the input filename stem (stripping known suffixes such as
        ``.cif``, ``.cif.gz``).

        Returns:
            :class:`pathlib.Path` to the written file.
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Derive a clean entry name from the input filename
        stem = self.cif_file.name
        for suffix in (".cif.gz", ".cif"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break

        out_path = self.output_dir / f"{stem}_edited.cif"
        self._doc.write_file(str(out_path))
        logger.info(f"Enriched CIF written to {out_path}")
        return out_path


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def run() -> None:
    """Command-line entry point: ``pdbe_sifts edit_cif``.

    Args (parsed from command line):
        -i / --input-cif:           Path to the input mmCIF file.
        -e / --entity-id:           Entity ID (e.g. ``1``).  Optional — when
                                    omitted (together with ``-c``) all
                                    polypeptide entity/chain pairs are
                                    processed automatically.
        -c / --chain-id:            Author chain ID (e.g. ``A``).  Optional —
                                    must be supplied together with ``-e`` or
                                    not at all.
        -o / --output-dir:          Output directory.
        --add-poly-seq-scheme:      Inject ``_pdbx_poly_seq_scheme`` (default:
                                    injected when absent).
        --no-poly-seq-scheme:       Skip poly_seq_scheme reconstruction.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="pdbe_sifts edit_cif",
        description="Enrich a CIF file by injecting missing mmCIF categories.",
    )
    parser.add_argument(
        "-i", "--input-cif", required=True, help="Input mmCIF file."
    )
    parser.add_argument(
        "-e",
        "--entity-id",
        default=None,
        help="Entity ID (e.g. 1). Omit to process all polypeptide entities.",
    )
    parser.add_argument(
        "-c",
        "--chain-id",
        default=None,
        help="Author chain ID (e.g. A). Omit to process all chains.",
    )
    parser.add_argument(
        "-o", "--output-dir", required=True, help="Output directory."
    )
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument(
        "--add-poly-seq-scheme",
        action="store_true",
        default=False,
        help="Force injection of _pdbx_poly_seq_scheme (even if already present).",
    )
    grp.add_argument(
        "--no-poly-seq-scheme",
        action="store_true",
        default=False,
        help="Skip _pdbx_poly_seq_scheme reconstruction entirely.",
    )
    args = parser.parse_args()

    # Validate: -e and -c must be supplied together or not at all
    if (args.entity_id is None) != (args.chain_id is None):
        parser.error(
            "Arguments -e/--entity-id and -c/--chain-id must be supplied "
            "together or both omitted."
        )

    ec = EditCif(
        args.input_cif,
        args.output_dir,
        entity_id=args.entity_id,
        chain_id=args.chain_id,
    )

    if not args.no_poly_seq_scheme:
        # Remove existing category when --add-poly-seq-scheme is forced
        if args.add_poly_seq_scheme:
            ec._block.find_mmcif_category(
                "_pdbx_poly_seq_scheme"
            )  # no-op if absent
        ec.add_pdbx_poly_seq_scheme()

    out = ec.write()
    print(f"Written: {out}")
