# Quick Start

This page walks through the complete SIFTS pipeline in six steps using the `pdbe_sifts` CLI.

## Step 1 — Initialise your config

```bash
pdbe_sifts init
# → creates ~/.config/pdbe_sifts/config.yaml
# → downloads the NCBI taxonomy database (~70 MB, first run only)
```

Edit the config to set your paths (`base_dir`, `nobackup_dir`, `target_db` after building it, etc.).

## Step 2 — Build a reference database

```bash
pdbe_sifts build_db \
  -i uniprot_sprot.fasta \
  -o ./my_db \
  -t taxonomy_mapping.tsv   # TSV: sequence_id <tab> tax_id
```

## Step 3 — Run global mappings

```bash
# Single CIF entry
pdbe_sifts global_mappings -i 1abc.cif -o ./results -d ./my_db/target_db

# Batch (one mmCIF path per line)
pdbe_sifts global_mappings -i entries.txt -o ./results -d ./my_db/target_db --threads 8
```

Produces `hits.duckdb` and `hits.tsv` — a scored table of UniProt accession candidates per PDB entity.

## Step 4 — Generate SIFTS segments and residue mappings

```bash
# With DuckDB hits (from global_mappings step)
pdbe_sifts segments -i 1abc.cif.gz -o ./segments -d hits.duckdb

# Manual UniProt accession mapping (chain:accession)
pdbe_sifts segments -i 1abc.cif.gz -o ./segments -m "A:P00963,B:P00963"

# Custom FASTA mapping (headers: >{auth_asym_id}|{sequence_id})
pdbe_sifts segments -i 1abc.cif.gz -o ./segments -m custom_seqs.fasta
```

Produces per-entry gzip-compressed CSV files under `{output_dir}/`.

## Step 5 — Load segment data into DuckDB

```bash
pdbe_sifts db_load -i ./segments/ -d hits.duckdb
```

Bulk-loads the segment and residue CSVs produced in step 4 into the `sifts_xref_segment` and `sifts_xref_residue` tables of the DuckDB file.

## Step 6 — Annotate mmCIF files with SIFTS data

```bash
# Reading segment data from DuckDB (after step 5)
pdbe_sifts sifts2mmcif \
  -i 1abc.cif.gz \
  -o ./sifts_mmcif \
  -d hits.duckdb

# Or reading segment CSVs directly from the output directory (skip step 5)
pdbe_sifts sifts2mmcif \
  -i 1abc.cif.gz \
  -o ./sifts_mmcif \
  -d hits.duckdb \
  -s ./segments/
```

!!! tip "Next steps"
    - See [CLI Reference](cli.md) for every flag on every subcommand.
    - See [API Reference](api/index.md) to drive the pipeline from Python without the CLI.
