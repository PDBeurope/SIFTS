# Quick Start

This page walks through the complete SIFTS pipeline in six steps using the `pdbe_sifts` CLI.

## Step 1 — Initialise your config

```bash
pdbe_sifts init
# → creates ~/.config/pdbe_sifts/config.yaml
# → downloads the NCBI taxonomy database (~70 MB, first run only)
# → builds the UniProt–PDB cross-reference index and CCD mapping cache
```

Edit the config to set your paths (`base_dir`, `nobackup_dir`, `target_db` after building it, etc.).

## Step 2 — Build a reference database

If you already have a plain or gzip-compressed UniProtKB FASTA, create its
taxonomy mapping and then build the database:

```bash
pdbe_sifts create_tax_file \
  --input-fasta uniprot_sprot.fasta.gz \
  --output-tax-mapping taxonomy_mapping.tsv

pdbe_sifts build_db \
  -i uniprot_sprot.fasta.gz \
  -o ./my_db/target_db \
  -t taxonomy_mapping.tsv \
  --tool mmseqs \
  --threads 8
```

`create_tax_file` expects UniProtKB `sp|...` or `tr|...` headers containing
`OX=<taxid>`. With custom headers, provide your own headerless
`sequence_id<TAB>tax_id` TSV to `build_db`.

## Step 3 — Run global mappings

```bash
# Single CIF entry
pdbe_sifts sequence_match -i 1abc.cif -o ./results -d ./my_db/target_db

# Batch (one mmCIF path per line)
pdbe_sifts sequence_match -i entries.txt -o ./results -d ./my_db/target_db --threads 8
```

Produces `hits.duckdb` and `hits_<entry>.tsv` under
`./results/mmseqs_<entry>/`.

## Step 4 — Generate SIFTS segments and residue mappings

```bash
# With DuckDB hits (from sequence_match step)
pdbe_sifts segments \
  -i 1abc.cif.gz \
  -o ./segments \
  -d ./results/mmseqs_1abc/hits.duckdb

# Custom FASTA mapping (headers: >{entry_id}|{auth_asym_id}|{sequence_id})
pdbe_sifts segments -i 1abc.cif.gz -o ./segments -m custom_seqs.fasta
```

Produces flat gzip-compressed files such as `{entry}_seg.csv.gz` and
`{entry}_res.csv.gz` under `{output_dir}/`.

## Step 5 — Load segment data into DuckDB

```bash
pdbe_sifts db_load \
  -i ./segments/ \
  -d ./results/mmseqs_1abc/hits.duckdb
```

Bulk-loads the segment and residue CSVs produced in step 4 into the `sifts_xref_segment` and `sifts_xref_residue` tables of the DuckDB file.

## Step 6 — Annotate mmCIF files with SIFTS data

```bash
# Reading segment data from DuckDB (after step 5)
pdbe_sifts sifts2mmcif \
  -i 1abc.cif.gz \
  -o ./sifts_mmcif \
  -d ./results/mmseqs_1abc/hits.duckdb

# Or reading segment CSVs directly from the output directory (skip step 5)
pdbe_sifts sifts2mmcif \
  -i 1abc.cif.gz \
  -o ./sifts_mmcif \
  -s ./segments/
```

!!! tip "Next steps"
    - See [CLI Reference](cli.md) for every flag on every subcommand.
    - See [API Reference](api/index.md) to drive the pipeline from Python without the CLI.
