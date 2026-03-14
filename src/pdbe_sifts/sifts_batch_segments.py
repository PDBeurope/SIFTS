"""Batch multiprocess segment generation for SIFTS.

Processes multiple CIF files in parallel using multiprocessing.Pool.
Each worker keeps a long-lived SiftsAlign instance (ChemCompMapping + DuckDB
loaded once) and swaps the CIF file path between entries.

Internal thread parallelism (lalign36 alignments, residue-map generation) is
controlled via the SIFTS_N_PROC environment variable so that
total_cpus ≈ n_workers × n_inner_threads.
"""

import os
from multiprocessing import Pool

import gemmi
from funcy.debug import log_durations

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import get_allocated_cpus
from pdbe_sifts.sifts_segments_generation import SiftsAlign

_state: dict = {}


def _derive_entry_id(cif_file: str) -> str:
    """Read _entry.id from a CIF file."""
    block = gemmi.cif.read(cif_file).sole_block()
    return block.find_value("_entry.id").strip('"').lower()


def _init_worker(first_cif, out_dir, db_conn_str, nf90_mode, unp_mode,
                 connectivity_mode, n_inner):
    """Create one long-lived SiftsAlign per worker process.

    ChemCompMapping and DuckDB connection are loaded once and reused for every
    entry the worker handles.
    """
    os.environ["SIFTS_N_PROC"] = str(n_inner)
    _state["align"] = SiftsAlign(
        first_cif, out_dir, db_conn_str, nf90_mode, unp_mode, connectivity_mode
    )


def _process_one(cif_file: str):
    """Process a single CIF file.  Called by ``pool.map()``."""
    try:
        align = _state["align"]
        align.cif_file = str(cif_file)
        align.custom_sequences = {}
        entry_id = _derive_entry_id(str(cif_file))
        align.process_entry(entry_id)
        return None
    except Exception:
        logger.error(f"Failed for {cif_file}", exc_info=True)
        return cif_file


@log_durations(logger.info)
def run_batch(cif_files, out_dir, db_conn_str=None, nf90_mode=False,
              unp_mode=None, connectivity_mode=True, workers=1):
    """Process multiple CIF files in parallel.

    Parameters
    ----------
    cif_files : list[str]
        Paths to input CIF files.
    out_dir : str
        Output directory for CSV files.
    db_conn_str : str | None
        DuckDB file path (optional when *unp_mode* is provided).
    nf90_mode : bool
        Enable UniRef90 mode.
    unp_mode : str | None
        User-defined mapping (accession string or FASTA path).
    connectivity_mode : bool
        Enable connectivity checking.
    workers : int
        Number of parallel worker processes.

    Returns
    -------
    list[str]
        List of CIF file paths that failed processing.
    """
    n_outer = workers
    n_inner = max(1, get_allocated_cpus() // n_outer)

    logger.info(
        f"Processing {len(cif_files)} entries with {n_outer} workers, "
        f"{n_inner} inner threads each"
    )

    first_cif = str(cif_files[0])

    with Pool(
        processes=n_outer,
        initializer=_init_worker,
        initargs=(first_cif, out_dir, db_conn_str, nf90_mode,
                  unp_mode, connectivity_mode, n_inner),
    ) as pool:
        results = pool.map(_process_one, [str(f) for f in cif_files])

    failed = [f for f in results if f is not None]
    if failed:
        logger.warning(f"{len(failed)} entries failed: {failed}")
    else:
        logger.info("All entries processed successfully")
    return failed
