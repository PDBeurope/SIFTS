#!/usr/bin/env python3

import argparse
from pathlib import Path
from pdbe_sifts.mmcif import extract_entities


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cif", required=True, help="Path to mmCIF file")
    parser.add_argument("--out", required=True, help="Output FASTA file")
    args = parser.parse_args()

    cif_path = Path(args.cif)
    pdb_id = cif_path.stem.lower()

    entities = extract_entities(str(cif_path))

    with open(args.out, "w") as fh:
        for entity_id, (sequence, taxid) in entities.items():
            header = f">pdb|{pdb_id}-{entity_id}|OX={taxid}"
            fh.write(header + "\n")
            fh.write(sequence + "\n")


if __name__ == "__main__":
    main()
