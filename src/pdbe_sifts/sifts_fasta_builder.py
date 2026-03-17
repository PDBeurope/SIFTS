#!/usr/bin/env python3
"""
FASTA builder for SIFTS.

Converts one or more mmCIF files into a single query FASTA file.
Supports three input types:
  - .fasta / .fa / .faa  — passed through as-is
  - .cif / .cif.gz       — single mmCIF file
  - .txt                  — text file listing mmCIF paths (one per line)
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.mmcif import extract_entities


class FastaBuilder:
    """Build a query FASTA file from mmCIF input(s) or an existing FASTA.

    Args:
        input_path: Path to a .cif/.cif.gz file, a .fasta/.fa/.faa file,
                    or a .txt file listing mmCIF paths (one per line).
        out_dir: Directory where the generated FASTA will be written.
                 Ignored when the input is already a FASTA file.
        threads: Number of parallel workers for .txt list processing.
        batch_size: Number of CIF files per batch for .txt list processing.
    """

    _FASTA_SUFFIXES = (".fasta", ".fa", ".faa")
    _CIF_SUFFIXES = (".cif.gz", ".cif")

    def __init__(
        self,
        input_path: str | Path,
        out_dir: str | Path,
        threads: int = 1,
        batch_size: int = 100000,
    ):
        self.input_path = Path(input_path)
        self.out_dir = Path(out_dir)
        self.threads = threads
        self.batch_size = batch_size

    @staticmethod
    def entry_name_from(input_path: str | Path) -> str:
        """Derive an entry name by stripping known file suffixes.

        Examples:
            1vyc.cif      → 1vyc
            1vyc.cif.gz   → 1vyc
            query.fasta   → query
            entries.txt   → entries
        """
        input_path = Path(input_path)
        name = input_path.name
        for suffix in (".cif.gz", ".cif", ".fasta", ".fa", ".faa", ".txt"):
            if name.lower().endswith(suffix):
                return name[: -len(suffix)]
        return input_path.stem

    @property
    def entry_name(self) -> str:
        """Entry name derived from the input filename (all known suffixes stripped)."""
        return self.entry_name_from(self.input_path)

    @staticmethod
    def _process_single_cif(file_path: str | Path):
        """Extract sequence/taxonomy from one mmCIF file.

        Returns:
            ((entry_name, entity_seq_tax_dict), None) on success
            (None, warning_message) on failure
        """
        file_path = Path(file_path)
        if not file_path.exists():
            return None, f"Missing CIF file: {file_path}"

        entry_name = FastaBuilder.entry_name_from(file_path)
        try:
            entity_seq_tax = extract_entities(file_path)
            return (entry_name, entity_seq_tax), None
        except Exception as e:
            return None, f"Failed to process {file_path}: {e}"

    @staticmethod
    def _write_fasta(
        entity_seq_tax_dict: dict,
        entry_name: str,
        fasta_path: Path,
        mode: str = "w",
    ) -> Path:
        """Write entity sequences to a FASTA file.

        Args:
            entity_seq_tax_dict: Mapping of entity_id → (sequence, taxid).
            entry_name: Entry name used in FASTA headers.
            fasta_path: Destination file path.
            mode: 'w' to create/overwrite, 'a' to append.

        Returns:
            Path to the written FASTA file.
        """
        with open(fasta_path, mode) as fh:
            for entity, (seq, tax_id) in entity_seq_tax_dict.items():
                fh.write(f">pdb|{entry_name}-{entity}|OX={tax_id}\n{seq}\n")
        return fasta_path

    def build(self) -> Path:
        """Build the FASTA file and return its path.

        Returns:
            Path to the FASTA file (existing or newly generated).

        Raises:
            FileNotFoundError: If the input file does not exist.
            ValueError: If the input format is not supported.
            RuntimeError: If no sequences could be extracted from a .txt input.
        """
        if not self.input_path.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_path}")

        name = self.input_path.name.lower()

        # FASTA input: pass through unchanged
        if name.endswith(self._FASTA_SUFFIXES):
            logger.info(f"Using FASTA input directly: {self.input_path}")
            return self.input_path

        # Single mmCIF file
        if name.endswith(self._CIF_SUFFIXES):
            logger.info(f"Extracting sequences from {self.input_path}")
            entity_seq_tax = extract_entities(self.input_path)
            fasta_path = self.out_dir / f"{self.entry_name}.fasta"
            self._write_fasta(entity_seq_tax, self.entry_name, fasta_path)
            return fasta_path

        # Text file listing mmCIF paths
        if name.endswith(".txt"):
            return self._build_from_list()

        raise ValueError(
            f"Unsupported input format: '{self.input_path.suffix}'. "
            "Expected .cif, .cif.gz, .fasta, .fa, .faa, or .txt"
        )

    def _build_from_list(self) -> Path:
        """Process a .txt file listing mmCIF paths, one per line."""
        fasta_path = self.out_dir / f"{self.entry_name}.fasta"
        fasta_path.unlink(missing_ok=True)

        file_paths = [
            line.strip()
            for line in self.input_path.read_text().splitlines()
            if line.strip()
        ]
        logger.info(
            f"Processing {len(file_paths)} CIF files from {self.input_path}"
        )

        for i in range(0, len(file_paths), self.batch_size):
            batch = file_paths[i : i + self.batch_size]
            results = []

            with ProcessPoolExecutor(max_workers=self.threads) as pool:
                futures = {
                    pool.submit(self._process_single_cif, fp): fp
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
                self._write_fasta(entity_seq_tax, en, fasta_path, mode="a")

        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
            raise RuntimeError(
                f"No sequences were extracted from {self.input_path}. "
                "Check that the CIF paths in the file are valid and readable."
            )

        return fasta_path


def run():
    """Command-line entry point for building a query FASTA from mmCIF input(s).

    Accepts a single .cif/.cif.gz, an existing .fasta/.fa/.faa,
    or a .txt file listing mmCIF paths (one per line).
    """
    parser = argparse.ArgumentParser(
        description=(
            "Build a query FASTA from one or more mmCIF files. "
            "Accepts a single .cif/.cif.gz, an existing .fasta, "
            "or a .txt file listing mmCIF paths."
        )
    )
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        help="Input file (.cif, .cif.gz, .fasta, .fa, .faa, or .txt list of CIF paths).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory where the generated FASTA will be written.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of parallel workers for .txt list processing (default: 1).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100000,
        help="Number of CIF files per batch for .txt list processing (default: 100000).",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = FastaBuilder(
        input_path=args.input_file,
        out_dir=out_dir,
        threads=args.threads,
        batch_size=args.batch_size,
    ).build()

    logger.info(f"FASTA written to: {fasta_path}")


if __name__ == "__main__":
    run()
