#!/usr/bin/env python3

import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from pdbe_sifts.mmcif import extract_entities


def process_cif(cif_path):
    cif_path = Path(cif_path)
    pdb_id = cif_path.stem.lower()

    entities = extract_entities(str(cif_path))

    records = []
    for entity_id, (sequence, taxid) in entities.items():
        header = f">pdb|{pdb_id}-{entity_id}|OX={taxid}"
        records.append(f"{header}\n{sequence}\n")

    return "".join(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--list", required=True, help="Text file with mmCIF paths (one per line)")
    parser.add_argument("--out", required=True, help="Merged FASTA output file")
    parser.add_argument(
        "--workers",
        type=int,
        default=cpu_count(),
        help="Number of processes (default: all CPUs)",
    )
    args = parser.parse_args()

    with open(args.list) as fh:
        cif_files = [line.strip() for line in fh if line.strip()]

    with Pool(processes=args.workers) as pool:
        results = pool.map(process_cif, cif_files)

    with open(args.out, "w") as out_fh:
        for fasta_chunk in results:
            out_fh.write(fasta_chunk)


if __name__ == "__main__":
    main()
