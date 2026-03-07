#!/usr/bin/env python3
"""
SIFTS global mappings generator.

This module performs global sequence–structure mappings for PDB entries
using MMseqs2 or BLASTP. It extracts sequences from mmCIF files, builds
FASTA files, runs alignments, and parses the resulting hits into ranked
mappings stored in DuckDB.
"""

import argparse
import shutil
from timeit import default_timer as timer
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

from pdbe_sifts.base.log import logger
from pdbe_sifts.mmcif import extract_entities
from pdbe_sifts.global_mappings.mmseqs_search import MmSearch
from pdbe_sifts.global_mappings.blastp import BlastP
from pdbe_sifts.global_mappings.global_mappings_parser import GlobMappingsParser
from pdbe_sifts.base.utils import get_date


class SiftsGlobalMappings:
    """Main class for generating structure–sequence mappings."""

    def __init__(
        self,
        input_file: str | Path,
        out_dir: str | Path,
        db_file: str | Path,
        unp_csv: Optional[str | Path] = None,
        tool: str = "mmseqs",
        threads: int = 1,
        batch_size: int = 100000,
    ):
        """
        Initialize SiftsGlobalMappings.

        Args:
            input_file: Path to a mmCIF file (.cif / .cif.gz), a FASTA file
                        (.fasta / .fa / .faa), or a text file listing mmCIF
                        paths (.txt).
            out_dir: Directory where results will be written.
            db_file: Path to the preformatted sequence database (MMseqs or BLAST).
            unp_csv: Optional CSV with accession metadata
                     (accession, provenance, pdb_xref, annotation_score).
            tool: Alignment tool ('mmseqs' or 'blastp'). Default: 'mmseqs'.
            threads: Number of CPU threads to use.
            batch_size: Number of CIF files per batch when processing a .txt list.
        """
        self.input_file = Path(input_file)
        self.out_dir = Path(out_dir)
        self.db_file = Path(db_file)
        self.unp_csv = Path(unp_csv) if unp_csv is not None else None
        self.tool = tool
        self.threads = threads
        self.batch_size = batch_size
        self.date = get_date()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def entry_name_from(input_file: str | Path) -> str:
        """Derive an entry name by stripping all known file suffixes.

        Examples:
            1vyc.cif      → 1vyc
            1vyc.cif.gz   → 1vyc
            query.fasta   → query
            entries.txt   → entries
        """
        input_file = Path(input_file)
        name = input_file.name
        for suffix in (".cif.gz", ".cif", ".fasta", ".fa", ".faa", ".txt"):
            if name.lower().endswith(suffix):
                return name[: -len(suffix)]
        return input_file.stem

    @property
    def entry_name(self) -> str:
        """Entry/run name derived from the input filename."""
        return self.entry_name_from(self.input_file)

    # ------------------------------------------------------------------
    # FASTA generation
    # ------------------------------------------------------------------

    @staticmethod
    def generate_fasta(
        entity_seq_tax_dict: dict,
        entry_name: str,
        fasta_path: Path,
        mode: str = "w",
    ) -> Path:
        """Write entity sequences to a FASTA file.

        Args:
            entity_seq_tax_dict: Mapping of entity_id → (sequence, taxid).
            entry_name: Entry name used in the FASTA header (e.g. '1vyc').
            fasta_path: Destination file path.
            mode: File open mode ('w' to create/overwrite, 'a' to append).

        Returns:
            Path to the written FASTA file.
        """
        with open(fasta_path, mode) as fh:
            for entity, (seq, tax_id) in entity_seq_tax_dict.items():
                fh.write(f">pdb|{entry_name}-{entity}|OX={tax_id}\n{seq}\n")
        return fasta_path

    # ------------------------------------------------------------------
    # Input processing
    # ------------------------------------------------------------------

    @staticmethod
    def _process_single_file(file_path: str | Path):
        """Extract sequence/taxonomy info from one mmCIF file.

        Module-level staticmethod so it can be pickled by ProcessPoolExecutor.

        Returns:
            ((entry_name, entity_seq_tax_dict), None) on success
            (None, warning_message) on failure
        """
        file_path = Path(file_path)
        if not file_path.exists():
            return None, f"Missing CIF file: {file_path}"

        name = file_path.name
        entry_name = name
        for suffix in (".cif.gz", ".cif"):
            if name.lower().endswith(suffix):
                entry_name = name[: -len(suffix)]
                break
        else:
            entry_name = file_path.stem

        try:
            entity_seq_tax = extract_entities(file_path)
            return (entry_name, entity_seq_tax), None
        except Exception as e:
            return None, f"Failed to process {file_path}: {e}"

    def process_input_file(self, result_dir: Path) -> tuple[Path, str]:
        """Process the input file and return (fasta_path, entry_name).

        Supports three input formats:
        - .fasta / .fa / .faa  — used as-is (no conversion)
        - .cif / .cif.gz       — single mmCIF, sequences extracted via extract_entities()
        - .txt                  — list of mmCIF paths, processed in parallel batches

        The generated FASTA (for CIF inputs) is written into result_dir so it
        lives alongside hits.duckdb.

        Args:
            result_dir: The run's result directory (already created by process()).

        Returns:
            Tuple of (path to FASTA file, entry/run name).
        """
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")

        name = self.input_file.name.lower()
        entry_name = self.entry_name

        # ---- FASTA input: use directly, no conversion needed ----
        if name.endswith((".fasta", ".fa", ".faa")):
            logger.info(f"Using FASTA input directly: {self.input_file}")
            return self.input_file, entry_name

        # ---- Single mmCIF file ----
        if name.endswith((".cif.gz", ".cif")):
            logger.info(f"Extracting sequences from {self.input_file}")
            entity_seq_tax = extract_entities(self.input_file)
            fasta_path = result_dir / f"{entry_name}.fasta"
            self.generate_fasta(entity_seq_tax, entry_name, fasta_path)
            return fasta_path, entry_name

        # ---- Text file listing mmCIF paths ----
        if name.endswith(".txt"):
            fasta_path = result_dir / f"{entry_name}.fasta"
            fasta_path.unlink(missing_ok=True)

            file_paths = [
                line.strip()
                for line in self.input_file.read_text().splitlines()
                if line.strip()
            ]
            logger.info(f"Processing {len(file_paths)} CIF files from {self.input_file}")

            for i in range(0, len(file_paths), self.batch_size):
                batch = file_paths[i : i + self.batch_size]
                results = []

                with ProcessPoolExecutor(max_workers=self.threads) as pool:
                    futures = {
                        pool.submit(self._process_single_file, fp): fp
                        for fp in batch
                    }
                    for fut in as_completed(futures):
                        res, warn = fut.result()
                        if warn:
                            logger.warning(warn)
                            continue
                        if res:
                            results.append(res)

                for en, entity_seq_tax in results:
                    self.generate_fasta(entity_seq_tax, en, fasta_path, mode="a")

            if not fasta_path.exists() or fasta_path.stat().st_size == 0:
                raise RuntimeError(
                    f"No sequences were extracted from {self.input_file}. "
                    "Check that the CIF paths in the file are valid and readable."
                )

            return fasta_path, entry_name

        raise ValueError(
            f"Unsupported input format '{self.input_file.suffix}'. "
            "Expected .cif, .cif.gz, .fasta, .fa, .faa, or .txt"
        )

    # ------------------------------------------------------------------
    # Search
    # ------------------------------------------------------------------

    def mmseqs_search(self, entry_id: str, fasta_path: Path, result_dir: Path) -> Path:
        """Run an MMseqs2 search."""
        output_path = result_dir / f"hits_{entry_id}.tsv"
        tmp_dir = result_dir / f"tmp_{entry_id}"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        MmSearch(
            fasta_path, self.db_file, output_path, tmp_dir,
            self.threads, db_load_mode=2,
        ).run()
        return output_path

    def blastp_search(self, entry_id: str, fasta_path: Path, result_dir: Path) -> Path:
        """Run a BLASTP search."""
        output_path = result_dir / f"hits_{entry_id}.tsv"
        BlastP(fasta_path, self.db_file, output_path, threads=self.threads).run()
        return output_path

    def search(self, entry_id: str, fasta_path: Path, result_dir: Path) -> Path:
        """Dispatch to the selected search tool."""
        if self.tool == "mmseqs":
            return self.mmseqs_search(entry_id, fasta_path, result_dir)
        elif self.tool == "blastp":
            return self.blastp_search(entry_id, fasta_path, result_dir)
        raise ValueError(f"Unsupported tool: {self.tool!r} (expected 'mmseqs' or 'blastp')")

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def process(self):
        """Run the complete mapping pipeline.

        Creates (or recreates) a result directory named
        ``{out_dir}/{tool}_{entry_name}_{date}/`` and writes all outputs
        there: the generated FASTA, the raw hits TSV, and hits.duckdb.
        """
        logger.info(f"Processing [{self.input_file}]")
        start = timer()

        # Create (or wipe and recreate) the result directory
        result_dir = self.out_dir / f"{self.tool}_{self.entry_name}_{self.date}"
        if result_dir.exists():
            shutil.rmtree(result_dir)
        result_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Results will be written to {result_dir}")

        fasta_path, entry_name = self.process_input_file(result_dir)
        hits_path = self.search(entry_name, fasta_path, result_dir)

        logger.info(f"Parsing {self.tool} hits.")
        GlobMappingsParser(
            self.tool,
            hits_path,
            result_dir,
            unp_csv=self.unp_csv,
        ).parse()

        end = timer()
        logger.info(f"Total time (extraction → ranked mappings): {end - start:.2f} s.")


def run():
    """Command-line entry point for generating SIFTS global mappings.

    Accepts a single mmCIF file, a FASTA file, or a text file listing
    mmCIF paths, runs an alignment against a reference sequence database
    (MMseqs2 or BLASTP), and writes ranked mappings to the output directory.
    """
    parser = argparse.ArgumentParser(
        description="Generate SIFTS mappings between a structure and sequence database."
    )
    parser.add_argument(
        "-i", "--input-file", required=True,
        help=(
            "Input file. Accepted formats: "
            "a single mmCIF file (.cif or .cif.gz), "
            "a FASTA file (.fasta / .fa / .faa), "
            "or a text file listing multiple mmCIF paths (.txt)."
        ),
    )
    parser.add_argument(
        "-o", "--output-dir", required=True,
        help="Directory where all results will be written.",
    )
    parser.add_argument(
        "-d", "--db-file", required=True,
        help="Path to the preformatted sequence database (MMseqs or BLAST).",
    )
    parser.add_argument(
        "--tool", default="mmseqs",
        help="Alignment tool to use ('mmseqs' or 'blastp'). Default: mmseqs.",
    )
    parser.add_argument(
        "--unp-csv-file", default=None,
        help=(
            "Path to CSV with accession metadata "
            "(accession, provenance, pdb_xref, annotation_score)."
        ),
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of threads to use for parsing and searches.",
    )
    parser.add_argument(
        "--batch-size", type=int, default=100000,
        help="Number of CIF files to process per batch when using a .txt list (default: 100000).",
    )

    args = parser.parse_args()
    logger.info(vars(args))

    pipeline = SiftsGlobalMappings(
        input_file=args.input_file,
        out_dir=args.output_dir,
        db_file=args.db_file,
        unp_csv=args.unp_csv_file,
        tool=args.tool,
        threads=args.threads,
        batch_size=args.batch_size,
    )
    pipeline.process()


if __name__ == "__main__":
    run()
