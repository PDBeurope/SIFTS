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
from pathlib import Path
from timeit import default_timer as timer

from pdbe_sifts.base.log import logger
from pdbe_sifts.sequence_match.blastp import BlastP
from pdbe_sifts.sequence_match.mmseqs_search import MmSearch
from pdbe_sifts.sequence_match.sequence_match_parser import SequenceMatchParser
from pdbe_sifts.sifts_fasta_builder import FastaBuilder


class SiftsSequenceMatch:
    """Main class for generating structure–sequence mappings."""

    def __init__(
        self,
        input_file: str | Path,
        out_dir: str | Path,
        db_file: str | Path,
        unp_csv: str | Path | None = None,
        tool: str = "mmseqs",
        threads: int = 1,
        batch_size: int = 100000,
    ):
        """
        Initialize SiftsSequenceMatch.

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

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @property
    def entry_name(self) -> str:
        """Entry/run name derived from the input filename (all known suffixes stripped)."""
        return FastaBuilder.entry_name_from(self.input_file)

    # ------------------------------------------------------------------
    # Input processing — delegated to FastaBuilder
    # ------------------------------------------------------------------

    def process_input_file(self, result_dir: Path) -> tuple[Path, str]:
        """Process the input file and return (fasta_path, entry_name).

        Delegates all FASTA generation logic to FastaBuilder, which supports:
        - .fasta / .fa / .faa  — used as-is (passthrough)
        - .cif / .cif.gz       — single mmCIF, sequences extracted
        - .txt                  — list of mmCIF paths, processed in parallel batches

        Args:
            result_dir: The run's result directory (already created by process()).

        Returns:
            Tuple of (path to FASTA file, entry/run name).
        """
        fasta_path = FastaBuilder(
            self.input_file, result_dir, self.threads, self.batch_size
        ).build()
        return fasta_path, FastaBuilder.entry_name_from(self.input_file)

    # ------------------------------------------------------------------
    # Search
    # ------------------------------------------------------------------

    def mmseqs_search(
        self, entry_id: str, fasta_path: Path, result_dir: Path
    ) -> Path:
        """Run an MMseqs2 easy_search against the configured database.

        Args:
            entry_id: Entry/run name used to name output files.
            fasta_path: Path to the query FASTA file.
            result_dir: Directory where TSV hits and temp files will be written.

        Returns:
            Path to the resulting hits TSV file.
        """
        output_path = result_dir / f"hits_{entry_id}.tsv"
        tmp_dir = result_dir / f"tmp_{entry_id}"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        MmSearch(
            fasta_path, self.db_file, output_path, tmp_dir, self.threads
        ).run()
        return output_path

    def blastp_search(
        self, entry_id: str, fasta_path: Path, result_dir: Path
    ) -> Path:
        """Run a BLASTP search against the configured BLAST database.

        Args:
            entry_id: Entry/run name used to name output files.
            fasta_path: Path to the query FASTA file.
            result_dir: Directory where the TSV hits file will be written.

        Returns:
            Path to the resulting hits TSV file.
        """
        output_path = result_dir / f"hits_{entry_id}.tsv"
        BlastP(
            fasta_path, self.db_file, output_path, threads=self.threads
        ).run()
        return output_path

    def search(self, entry_id: str, fasta_path: Path, result_dir: Path) -> Path:
        """Dispatch to the search tool selected at construction time.

        Args:
            entry_id: Entry/run name used to name output files.
            fasta_path: Path to the query FASTA file.
            result_dir: Directory where output files will be written.

        Returns:
            Path to the resulting hits TSV file.

        Raises:
            ValueError: If ``self.tool`` is not ``"mmseqs"`` or ``"blastp"``.
        """
        if self.tool == "mmseqs":
            return self.mmseqs_search(entry_id, fasta_path, result_dir)
        elif self.tool == "blastp":
            return self.blastp_search(entry_id, fasta_path, result_dir)
        raise ValueError(
            f"Unsupported tool: {self.tool!r} (expected 'mmseqs' or 'blastp')"
        )

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def process(self) -> None:
        """Run the complete mapping pipeline.

        Creates (or recreates) a result directory named
        ``{out_dir}/{tool}_{entry_name}_{date}/`` and writes all outputs
        there: the generated FASTA, the raw hits TSV, and hits.duckdb.
        """
        logger.info(f"Processing [{self.input_file}]")
        start = timer()

        # Create (or wipe and recreate) the result directory
        result_dir = self.out_dir / f"{self.tool}_{self.entry_name}"
        if result_dir.exists():
            shutil.rmtree(result_dir)
        result_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Results will be written to {result_dir}")

        fasta_path, entry_name = self.process_input_file(result_dir)
        hits_path = self.search(entry_name, fasta_path, result_dir)

        logger.info(f"Parsing {self.tool} hits.")
        SequenceMatchParser(
            self.tool,
            hits_path,
            result_dir,
            unp_csv=self.unp_csv,
        ).parse()

        end = timer()
        logger.info(
            f"Total time (extraction → ranked mappings): {end - start:.2f} s."
        )


def run() -> None:
    """Command-line entry point for generating SIFTS global mappings.

    Accepts a single mmCIF file, a FASTA file, or a text file listing
    mmCIF paths, runs an alignment against a reference sequence database
    (MMseqs2 or BLASTP), and writes ranked mappings to the output directory.
    """
    parser = argparse.ArgumentParser(
        description="Generate SIFTS mappings between a structure and sequence database."
    )
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        help=(
            "Input file. Accepted formats: "
            "a single mmCIF file (.cif or .cif.gz), "
            "a FASTA file (.fasta / .fa / .faa), "
            "or a text file listing multiple mmCIF paths (.txt)."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory where all results will be written.",
    )
    parser.add_argument(
        "-d",
        "--db-file",
        required=True,
        help="Path to the preformatted sequence database (MMseqs or BLAST).",
    )
    parser.add_argument(
        "--tool",
        default="mmseqs",
        help="Alignment tool to use ('mmseqs' or 'blastp'). Default: mmseqs.",
    )
    parser.add_argument(
        "--unp-csv-file",
        default=None,
        help=(
            "Path to CSV with accession metadata "
            "(accession, provenance, pdb_xref, annotation_score)."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for parsing and searches.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100000,
        help="Number of CIF files to process per batch when using a .txt list (default: 100000).",
    )

    args = parser.parse_args()
    logger.info(vars(args))

    pipeline = SiftsSequenceMatch(
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
