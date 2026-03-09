#!/usr/bin/env python3
"""Entry point executed on each LSF/Slurm compute node.

Usage::

    python process_remote.py /path/to/pickle_file

The script unpickles the :class:`~pdbe_sifts.base.batchable.Batchable`
subclass instance and calls its :meth:`~pdbe_sifts.base.batchable.Batchable.process_remote`
method, which pops entries from the shared queue and processes them.
"""

import pickle  # nosec
import sys


def run() -> None:
    """Main entry point — load pickle and call process_remote()."""
    if len(sys.argv) != 2:
        raise SystemExit(
            f"Usage: {sys.argv[0]} </path/to/pickle_file>\n"
            f"Got: {sys.argv}"
        )
    obj = _load_pickle(sys.argv[1])
    obj.process_remote()


def _load_pickle(file_name: str) -> object:
    """Deserialise and return the object stored in *file_name*."""
    with open(file_name, "rb") as fh:
        return pickle.load(fh)  # nosec


if __name__ == "__main__":
    run()
