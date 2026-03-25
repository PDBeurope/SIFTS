#!/usr/bin/env python3
import argparse
import gzip
import shutil
import tempfile
from pathlib import Path

import duckdb
from gemmi import cif

from pdbe_sifts.base.exceptions import EntryFailedException, NoSegmentsError
from pdbe_sifts.base.log import logger
from pdbe_sifts.sifts_to_mmcif import comm_utils
from pdbe_sifts.sifts_to_mmcif.def_mmcif_cat import (
    NEW_MMCIF_CAT,
    PRIMARY_KEYS,
    SIFTS_ATOMSITE_ITEM,
    SIFTS_NEW_CAT,
)
from pdbe_sifts.sifts_to_mmcif.delta_mappings import (
    FindMappingChanges,
    get_delta_csv_suffix,
    read_sifts_segments,
)
from pdbe_sifts.sifts_to_mmcif.read_sifts_csv import (
    get_unp_segments,
    get_unpres_mapping,
)


class ExportSIFTSTommCIF:
    """Write the updated mmcif file with sifts_data."""

    def __init__(
        self,
        input_cif,
        output_dir,
        sifts_csv_dir,
        duckdb_file,
        track_changes,
        prev_run_dir,
    ):
        """Initialise the SIFTS-to-mmCIF exporter.

        Args:
            input_cif (str | Path): Path to the input ``*_updated.cif.gz`` file.
            output_dir (str | Path): Directory where output mmCIF files will be
                written.
            sifts_csv_dir (str | Path | None): Directory containing pre-computed
                ``{entry_id}_seg.csv.gz`` and ``{entry_id}_res.csv.gz`` files.
                Pass ``None`` to read data from the DuckDB database instead.
            duckdb_file (str | Path | None): Path to the DuckDB database file.
                Pass ``None`` when ``sifts_csv_dir`` is provided for all entries.
            track_changes (bool): When ``True``, compare the new SIFTS output
                against the previous run and write a delta CSV file.
            prev_run_dir (str | Path | None): Directory containing
                ``{entry_id}_sifts_only.cif.gz`` from the previous run, used
                as the baseline for change-tracking.  When ``None`` the output
                directory is used as the baseline.
        """
        self.input_cif = input_cif
        self.output_dir = output_dir
        self.duckdb_file = duckdb_file
        self.sifts_csv_base = sifts_csv_dir
        self.conn = (
            duckdb.connect(self.duckdb_file, read_only=True)
            if duckdb_file
            else None
        )
        self.track_changes = track_changes
        self.prev_run_dir = prev_run_dir

    def _check_mmcif_output(self, infile, outfile):
        """Verify that the output mmCIF file is at least as large as the input.

        If the output is smaller than the input, a deeper content-level check
        is performed via :func:`comm_utils.check_output_mmcif`.

        Args:
            infile (Path): Path to the original input CIF file.
            outfile (Path): Path to the generated output CIF file.

        Raises:
            EntryFailedException: If ``outfile`` does not exist.
        """
        logger.debug("Checking if the output file is more/equal to input file")

        if outfile.exists():
            input_size = infile.stat().st_size
            output_size = outfile.stat().st_size
            if output_size < input_size:
                logger.warning(
                    "Check output file- it's smaller than input) for "
                    f"{outfile} opsize vs insize --> {output_size} vs {input_size}"
                )
                comm_utils.check_output_mmcif(infile, outfile)

        else:
            raise EntryFailedException(
                f"Output {outfile} file not found for {outfile}"
            )

    def _check_clean_mmcif(self, outfile):
        """Check that the output clean mmcif file has no null values in primary keys."""
        block = cif.read(str(outfile)).sole_block()
        category_list = [f"_{item}" for item in SIFTS_NEW_CAT]
        required_fields = PRIMARY_KEYS
        required_fields.setdefault("_atom_site", []).append(
            [
                "pdbx_sifts_xref_db_name",
                "pdbx_sifts_xref_db_acc",
                "pdbx_sifts_xref_db_num",
                "pdbx_sifts_xref_db_res",
            ]
        )
        for cat in category_list:
            cols = required_fields[cat]
            for col in cols:
                if any(
                    cif.is_null(x) for x in block.find_values(f"{cat}.{col}")
                ):
                    raise Exception(
                        f"Data error invalid value (? or .) for required item: {cat}.{col}"
                    )

    def _check_sifts_mmcif(self, sifts_only_cif: Path):
        """Verify that the SIFTS-only mmCIF output file is non-empty.

        Args:
            sifts_only_cif (Path): Path to the ``*_sifts_only.cif.gz`` file.

        Raises:
            EntryFailedException: If the file does not exist or its size is
                40 bytes or fewer (considered effectively empty).
        """
        logger.debug("Checking sifts mmcif output is not empty")
        if sifts_only_cif.exists():
            if sifts_only_cif.stat().st_size <= 40:
                raise EntryFailedException(
                    f"Empty sifts output file {sifts_only_cif}"
                )
        else:
            raise EntryFailedException(
                f"file not found sifts mmcif file {sifts_only_cif}"
            )

    def _process_cif_file(
        self, d_block, sifts_segments, sifts_seg_inst, sifts_res_csv, pdb_id
    ):
        """Populate a CIF block with SIFTS segment and residue-level data.

        Writes the ``_pdbx_sifts_unp_segments`` and ``_pdbx_sifts_xref_db``
        categories into ``d_block`` and updates ``_atom_site`` with the SIFTS
        cross-reference columns.

        Args:
            d_block: A ``gemmi`` CIF block to be updated in place.
            sifts_segments: Column data (tuple of lists) for
                ``_pdbx_sifts_unp_segments``, as returned by
                :func:`read_sifts_csv.get_unp_segments`.
            sifts_seg_inst (dict): Segment-instance lookup used to expand
                segment mappings to individual residues.
            sifts_res_csv (Path | None): Path to the residue-level CSV file.
                ``None`` when data is read from the DuckDB connection.
            pdb_id (str): The PDB entry identifier.

        Returns:
            tuple: ``(UpdateFlag, d_block, UnpDataFlag)`` where ``UpdateFlag``
            is ``True`` if any SIFTS data was written, ``d_block`` is the
            modified CIF block, and ``UnpDataFlag`` is ``True`` when UniProt
            residue data was successfully applied to ``_atom_site``.
        """
        self.UpdateFlag = False
        self.UnpDataFlag = False

        sifts_res_info = {}

        self.UpdateFlag = self.write_data(
            d_block, "_pdbx_sifts_unp_segments", sifts_segments
        )

        cat = d_block.get_mmcif_category("_pdbx_poly_seq_scheme")
        if cat:
            my_ent = cat["entity_id"]
            my_ch = cat["asym_id"]
            my_seq_id = cat["seq_id"]
            my_mon_id = cat["mon_id"]
            my_observed = cat["auth_seq_num"]
            my_observed = [
                "y" if item is not None else "n" for item in my_observed
            ]
            my_obs_res = comm_utils.get_obs(my_ch, my_seq_id, my_observed)
            hetero = cat["hetero"]
            my_mh_id = comm_utils.get_mh_id(my_seq_id, my_mon_id, hetero)

            all_ent, all_res = comm_utils.get_ent_chains(
                my_ent, my_ch, my_seq_id, my_mon_id
            )

            if sifts_res_csv:
                res_cursor = None
            else:
                res_cursor = self.conn
            sifts_res_info, mon_id_info = get_unpres_mapping(
                pdb_id, sifts_res_csv, res_cursor
            )
            xref_resmap = comm_utils.get_xref_db(
                all_res,
                sifts_res_info,
                sifts_seg_inst,
                my_obs_res,
                my_mh_id,
                mon_id_info,
            )
            xref_resmap = [item for item in xref_resmap if item != []]

            self.UpdateFlag = self.write_data(
                d_block, "_pdbx_sifts_xref_db", xref_resmap
            )

        else:
            logger.warning(f"The entry {pdb_id} is not a polymer")

        if sifts_res_info:
            self.UnpDataFlag = True
            self.UpdateFlag = True
            atom_site = d_block.get_mmcif_category("_atom_site")
            atom_site = {
                "pdbx_label_index" if k == "pdbe_label_seq_id" else k: v
                for k, v in atom_site.items()
            }
            atom_site = comm_utils.modify_atomsite(atom_site, sifts_res_info)
            d_block.set_mmcif_category("_atom_site", atom_site)
        else:
            if not self.UnpDataFlag:
                atom_site = d_block.get_mmcif_category("_atom_site")
                atom_site = {
                    "pdbx_label_index" if k == "pdbe_label_seq_id" else k: v
                    for k, v in atom_site.items()
                }
                d_block.set_mmcif_category("_atom_site", atom_site)

        return self.UpdateFlag, d_block, self.UnpDataFlag

    def write_data(self, d_block, category, category_data):
        """Write a SIFTS mmCIF category into a CIF block.

        The data is written only when ``category_data`` is non-empty and its
        length matches the number of columns defined for ``category`` in
        :data:`NEW_MMCIF_CAT`.

        Args:
            d_block: The ``gemmi`` CIF block to update.
            category (str): The mmCIF category name (e.g.
                ``"_pdbx_sifts_unp_segments"``).
            category_data (list | None): Parallel column lists to be written.
                Must have the same length as ``NEW_MMCIF_CAT[category]``.

        Returns:
            bool: ``True`` if the category was successfully written, ``False``
            otherwise.
        """
        if category_data and len(NEW_MMCIF_CAT[category]) == len(category_data):
            data = dict(
                zip(NEW_MMCIF_CAT[category], category_data, strict=False)
            )
            d_block.set_mmcif_category(category, data)
            return True
        return False

    def delta_changes(self):
        """Compare new SIFTS output against the previous run and record changes.

        When change-tracking is enabled (``self.track_changes``) and UniProt
        data was written (``self.UnpDataFlag``), this method determines the
        delta CSV suffix, reads the old and new segment data, and invokes
        :class:`~pdbe_sifts.sifts_to_mmcif.delta_mappings.FindMappingChanges`
        to write the delta CSV file.  If no changes are detected, no output
        file is produced.

        Sets ``self.sifts_delta_csv`` to the output path when a delta file is
        written, or to ``None`` when tracking is skipped or no changes exist.
        """
        self.sifts_delta_csv = None
        if self.track_changes and self.UnpDataFlag:
            base_sifts_dir = Path(self.sifts_updated_cif).parent
            delta_suff = get_delta_csv_suffix(base_sifts_dir)
            self.sifts_delta_csv = Path(
                base_sifts_dir,
                f"{self.entry_id}_delta_{delta_suff}.csv.gz",
            )
            new_seg_data, new_seg_list = read_sifts_segments(self.d_block)
            if self.prev_run_dir:
                old_sifts_cif = (
                    Path(self.prev_run_dir)
                    / f"{self.entry_id}_sifts_only.cif.gz"
                )
            else:
                old_sifts_cif = self.sifts_only_cif

            old_seg_data, old_seg_list = {}, []
            if old_sifts_cif.exists():
                old_block = cif.read(str(old_sifts_cif)).sole_block()
                old_seg_data, old_seg_list = read_sifts_segments(old_block)
            else:
                logger.warning(
                    f"Track changes failed for {self.entry_id}."
                    " Previous run file not found!"
                    " Either a new entry or previous run failed!"
                )
            if old_seg_list == new_seg_list:
                logger.info(
                    f"No mappings changes tracked for {self.entry_id}."
                    " Skipping writing new mmcif files"
                )
                self.sifts_delta_csv = None
                return
            else:
                logger.info(f"Tracked mappings changed for {self.entry_id}")
                FindMappingChanges(
                    self.entry_id,
                    old_seg_data,
                    new_seg_data,
                    self.sifts_delta_csv,
                ).find_mapping_changes()

    def process_entry(self, entry_id):
        """Process a single PDB entry end-to-end.

        Reads the input CIF, merges SIFTS segment and residue data, writes the
        updated and SIFTS-only output files, runs change-tracking, and performs
        sanity checks on all outputs.

        Args:
            entry_id (str): The lowercase PDB entry identifier (e.g. ``"1abc"``).

        Returns:
            tuple[Path, Path, Path]: ``(orig_updated_cif, sifts_updated_cif,
            sifts_only_cif)`` — paths to the original input, the updated output,
            and the SIFTS-only output file.

        Raises:
            FileNotFoundError: If the input CIF file does not exist.
            NoSegmentsError: If no SIFTS segments could be written for the entry.
            EntryFailedException: If any output file fails its post-write checks.
        """
        out = Path(self.output_dir)
        out.mkdir(parents=True, exist_ok=True)
        self.sifts_updated_cif = out / f"{entry_id}_sifts_updated.cif.gz"
        self.sifts_only_cif = out / f"{entry_id}_sifts_only.cif.gz"

        if self.sifts_csv_base:
            seg_cursor = None
            sifts_seg_csv = Path(self.sifts_csv_base) / f"{entry_id}_seg.csv.gz"
            if not sifts_seg_csv.exists():
                logger.error(f"{sifts_seg_csv} does not exist")

            sifts_res_csv = Path(self.sifts_csv_base) / f"{entry_id}_res.csv.gz"
            if not sifts_res_csv.exists():
                logger.error(f"{sifts_res_csv} does not exist")
        else:
            sifts_seg_csv, sifts_res_csv = None, None
            seg_cursor = self.conn

        sifts_segments, sifts_seg_inst = get_unp_segments(
            entry_id, sifts_seg_csv, seg_cursor
        )
        sifts_segments = [item for item in sifts_segments if item != []]

        orig_updated_cif = Path(self.input_cif)

        if not orig_updated_cif.exists():
            raise FileNotFoundError(orig_updated_cif)

        self.d_block = cif.read(str(orig_updated_cif)).sole_block()

        self.UpdateFlag, self.d_block, self.UnpDataFlag = (
            self._process_cif_file(
                self.d_block,
                sifts_segments,
                sifts_seg_inst,
                sifts_res_csv,
                entry_id,
            )
        )

        self.entry_id = entry_id
        self.delta_changes()

        doc = cif.Document()
        doc.add_copied_block(self.d_block)
        with gzip.open(self.sifts_updated_cif, "wt") as f:
            f.write(doc.as_string(cif.Style.Pdbx))

        doc = cif.Document()
        block = doc.add_new_block(entry_id.upper())
        cat_to_copy = SIFTS_NEW_CAT

        for cat in cat_to_copy:
            data = self.d_block.get_mmcif_category(f"_{cat}")
            if not data:
                if cat == "_pdbx_sifts_unp_segments":
                    raise NoSegmentsError(
                        "No segments written. This shouldn't happen"
                    )
                logger.warning(f"Category _{cat} not written out")
                continue
            block.set_mmcif_category(f"_{cat}.", data)

        atom_cat = block.get_mmcif_category("_atom_site")
        if self.UnpDataFlag:
            for key in list(atom_cat.keys()):
                if key not in SIFTS_ATOMSITE_ITEM:
                    del atom_cat[key]
        else:
            for key in list(atom_cat.keys()):
                if key not in ["id", "pdbx_label_index"]:
                    del atom_cat[key]
        block.set_mmcif_category("_atom_site.", atom_cat)

        with tempfile.TemporaryDirectory() as d:
            temp_name = Path(d, f"{entry_id}_sifts_only.cif.gz")
            with gzip.open(temp_name, "wt") as f:
                f.write(doc.as_string(cif.Style.Pdbx))
                f.flush()
            shutil.copyfile(temp_name, self.sifts_only_cif)

        self._check_sifts_mmcif(self.sifts_only_cif)
        del self.d_block
        del block
        del doc

        if self.UpdateFlag:
            logger.info(f"SIFTS mapping found for {self.entry_id}")
        else:
            logger.info(f"No mappings in SIFTS for {self.entry_id}")

        self._check_mmcif_output(orig_updated_cif, self.sifts_updated_cif)
        self._check_clean_mmcif(self.sifts_updated_cif)
        logger.info(
            f"""
            Input CIF: {orig_updated_cif}
            Output files:
            UPDATED_MMCIF: {self.sifts_updated_cif}
            SIFTS_ONLY: {self.sifts_only_cif}
            DELTA_MAPPINGS: {self.sifts_delta_csv}
        """
        )
        return (orig_updated_cif, self.sifts_updated_cif, self.sifts_only_cif)


def run():
    """CLI entry point: parse arguments and run the SIFTS-to-mmCIF export.

    Parses command-line arguments, resolves the PDB entry identifier from the
    CIF file when not explicitly provided, constructs an
    :class:`ExportSIFTSTommCIF` instance, and calls
    :meth:`ExportSIFTSTommCIF.process_entry`.  The DuckDB connection is always
    closed in a ``finally`` block.

    Command-line arguments:
        --entry: Optional PDB entry ID.  Derived from ``_entry.id`` in the CIF
            when omitted.
        -i / --input-cif: Path to the input ``*_updated.cif.gz`` file.
        -o / --output-dir: Directory for output files.
        -s / --sifts-csv-dir: Optional directory containing CSV SIFTS data.
        -d / --db-conn-str: DuckDB file path.
        -T / --no-track-changes: Disable mapping-change tracking.
        -p / --prev-run-dir: Directory with the previous-run SIFTS-only CIF.
    """
    parser = argparse.ArgumentParser(
        "Given a clean mmcif file, add sifts_data",
    )
    parser.add_argument(
        "--entry",
        required=False,
        default=None,
        help="PDB entry ID to process. If omitted, derived from _entry.id in the CIF.",
    )
    parser.add_argument(
        "-i",
        "--input-cif",
        required=True,
        help="Input CIF file to process (*_updated.cif.gz)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Output directory where SIFTS mmcif files will be written",
    )
    parser.add_argument(
        "-s",
        "--sifts-csv-dir",
        required=False,
        help=(
            "Flat directory containing {entry_id}_seg.csv.gz / _res.csv.gz files. "
            "If not given, data is read from the database."
        ),
    )
    parser.add_argument(
        "-d",
        "--db-conn-str",
        required=True,
        help="DuckDB file path",
    )
    parser.add_argument(
        "-T",
        "--no-track-changes",
        required=False,
        action="store_true",
        default=False,
        help="Disable comparison with previous run to track mapping changes.",
    )
    parser.add_argument(
        "-p",
        "--prev-run-dir",
        required=False,
        help=(
            "Compare sifts_only.mmcif from this directory. "
            "If not given, looks in the output directory for the previous run file."
        ),
    )
    args = parser.parse_args()

    if not Path(args.input_cif).is_file():
        parser.error(f"-i must be a CIF file, got: {args.input_cif}")

    entry_id = args.entry
    if not entry_id:
        block = cif.read(str(args.input_cif)).sole_block()
        entry_id = block.find_value("_entry.id").strip('"').lower()

    track_changes = not args.no_track_changes
    obj = ExportSIFTSTommCIF(
        args.input_cif,
        args.output_dir,
        args.sifts_csv_dir,
        args.db_conn_str,
        track_changes,
        args.prev_run_dir,
    )
    try:
        obj.process_entry(entry_id)
    finally:
        obj.conn.close()


if __name__ == "__main__":
    run()
