# CLI Reference

All commands are invoked as:

```
pdbe_sifts [--log-level LEVEL] <command> [options]
```

## Global flag

| Flag | Default | Description |
|------|---------|-------------|
| `--log-level LEVEL` | `INFO` | Logging verbosity: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`. Overrides the `SIFTS_LOG_LEVEL` environment variable. |
| `--version` | — | Print the installed `pdbe_sifts` version and exit. |

## Command overview

| Command | Description |
|---------|-------------|
| [`init`](#init) | Create the config and initialise the taxonomy, UniProt–PDB, and CCD caches |
| [`show`](#show) | Print the fully resolved configuration |
| [`update_ncbi`](#update_ncbi) | Force-update the local NCBI taxonomy database (ete4) |
| [`prepare_build_db`](#prepare_build_db) | Prepare FASTA and taxonomy mapping inputs for database creation |
| [`create_tax_file`](#create_tax_file) | Extract an accession-to-taxid TSV from a UniProtKB FASTA |
| [`build_db`](#build_db) | Build a MMseqs2 or BLAST reference sequence database from a FASTA file |
| [`fasta_build`](#fasta_build) | Extract entity sequences from mmCIF files and write a FASTA |
| [`sequence_match`](#sequence_match) | Align PDB sequences against the reference DB and score hits into DuckDB |
| [`segments`](#segments) | Generate SIFTS segment and residue mappings for a **single** mmCIF entry |
| [`db_load`](#db_load) | Bulk-load segment/residue CSVs from the `segments` step into DuckDB |
| [`sifts2mmcif`](#sifts2mmcif) | Inject SIFTS mappings into an annotated mmCIF file |
| [`update_ccd_mapping`](#update_ccd_mapping) | Refresh the cached CCD three-to-one letter mapping when needed |
| [`seq2seq`](#seq2seq) | Align canonical deposited sequences against coordinate sequences |

---

## `init`

Copy the built-in config template to the user config directory, initialise the
NCBI taxonomy database, download the SIFTS UniProt–PDB cross-reference file,
build its local DuckDB index, and generate the CCD mapping cache.

```bash
pdbe_sifts init [--dest PATH] [--force] [--base-dir PATH] [--nobackup-dir PATH] [--target-db PATH]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--dest PATH` | `PDBE_SIFTS_CONFIG` or `~/.config/pdbe_sifts/config.yaml` | Custom destination path for the config file |
| `--force` | `False` | Overwrite an existing config file |
| `--base-dir PATH` | unset | Directory for PDBe-SIFTS result files |
| `--nobackup-dir PATH` | unset | Directory for large cache and temporary files |
| `--target-db PATH` | unset | Reference search database used by later commands |

---

## `show`

Print the fully resolved configuration (defaults merged with user overrides).

```bash
pdbe_sifts show [--config PATH]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--config PATH` | `~/.config/pdbe_sifts/config.yaml` | Path to a custom config file |

---

## `prepare_build_db`

Prepare a FASTA and its taxonomy mapping file, downloading UniProtKB/Swiss-Prot
when no input FASTA is supplied.

```bash
pdbe_sifts prepare_build_db --output-fasta FASTA.gz --output-tax-mapping TAX.tsv [--input-fasta FASTA] [--force]
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--output-fasta` | ✓ | — | Destination FASTA path |
| `--output-tax-mapping` | ✓ | — | Destination accession-to-taxid TSV path |
| `--input-fasta` | | downloaded Swiss-Prot | Existing plain or gzip-compressed UniProtKB FASTA |
| `--force` | | `False` | Overwrite an existing output FASTA |

## `create_tax_file`

Extract taxonomy mappings from an existing plain or gzip-compressed UniProtKB
FASTA without copying or downloading it.

```bash
pdbe_sifts create_tax_file --input-fasta FASTA[.gz] --output-tax-mapping TAX.tsv
```

| Flag | Required | Description |
|------|----------|-------------|
| `--input-fasta` | ✓ | Input UniProtKB `.fasta` or `.fasta.gz` file |
| `--output-tax-mapping` | ✓ | Output TSV path |

The TSV contains `ACCESSION<TAB>TAXID` rows without a header. Unsupported FASTA
headers are skipped; the command fails if no compatible headers are found.
Compatible headers begin with `sp|` or `tr|` and contain `OX=<taxid>`.

For an existing UniProtKB FASTA, create the taxonomy mapping and build the
database directly:

```bash
pdbe_sifts create_tax_file \
  --input-fasta FASTA.gz \
  --output-tax-mapping TAX.tsv

pdbe_sifts build_db \
  -i FASTA.gz \
  -o OUTPUT_PATH \
  -t TAX.tsv \
  --tool mmseqs \
  --threads 8
```

For custom FASTA headers, create the headerless `sequence_id<TAB>tax_id` TSV
yourself instead of using `create_tax_file`.

---

## `build_db`

Build an indexed reference sequence database from a FASTA file.

```bash
pdbe_sifts build_db -i FASTA -o OUTPUT_PATH -t TAX_FILE [--tool TOOL] [--threads N]
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i`, `--input-file` | ✓ | — | Input FASTA file |
| `-o`, `--output-path` | ✓ | — | Output path for database files |
| `-t`, `--tax-mapping-file` | ✓ | — | TSV file mapping sequence IDs to NCBI taxon IDs |
| `--tool` | | `mmseqs` | Search back-end: `mmseqs` or `blastp` |
| `--threads` | | `1` | Number of CPU threads |

MMseqs2 accepts gzip-compressed FASTA input. When using `--tool blastp`, provide
an input format supported by `makeblastdb`.

---

## `fasta_build`

Extract entity-level amino acid sequences from mmCIF files and write a multi-FASTA file.

```bash
pdbe_sifts fasta_build -i INPUT -o OUTPUT_DIR [--threads N] [--batch-size N]
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i`, `--input-file` | ✓ | — | `.cif`, `.cif.gz`, `.fasta`, `.fa`, `.faa`, or `.txt` listing CIF paths |
| `-o`, `--output-dir` | ✓ | — | Directory where the FASTA will be written |
| `--threads` | | `1` | Parallel workers for `.txt` list processing |
| `--batch-size` | | `100000` | CIF files per batch for `.txt` list processing |

---

## `sequence_match`

Align PDB entity sequences against the reference database, score hits with the SIFTS scoring function, and store results in DuckDB.

```bash
pdbe_sifts sequence_match -i INPUT [-o OUTPUT_DIR] [-d DB_FILE] [options]
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i`, `--input-file` | ✓ | — | `.cif`, `.cif.gz`, `.fasta`, `.fa`, `.faa`, or `.txt` listing CIF paths |
| `-o`, `--output-dir` | | config `base_dir` | Output directory for results |
| `-d`, `--db-file` | | config `target_db` | Pre-built reference database path |
| `--tool` | | `mmseqs` | Alignment tool: `mmseqs` or `blastp` |
| `--unp-csv-file` | | `None` | CSV with accession metadata (accession, provenance, annotation_score, …) |
| `--threads` | | `1` | Number of CPU threads |
| `--batch-size` | | `100000` | CIF files per batch for `.txt` list processing |
| `--tax-tsv` | | `None` | Query taxonomy overrides with columns `entry`, `entity`, `tax_id` |

**Outputs:** `hits.duckdb` and `hits_<entry>.tsv` under
`OUTPUT_DIR/{tool}_{entry_name}/` (for example, `OUTPUT_DIR/mmseqs_1abc/`).

---

## `segments`

Generate residue-level SIFTS segment and residue mappings for a single mmCIF entry.

```bash
pdbe_sifts segments -i CIF -o OUTPUT_DIR [-d DB] [-m MAPPING_FASTA] [options]
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i`, `--input-cif` | ✓ | — | Input CIF file (`.cif` or `.cif.gz`) |
| `-o`, `--output-dir` | ✓ | — | Output directory for per-entry CSV files |
| `-d`, `--db-file` | ✗* | — | DuckDB hits file from `sequence_match` |
| `-m`, `--mapping` | ✗* | — | Existing custom FASTA with headers `>{entry_id}|{auth_asym_id}|{sequence_id}[|{name}]` |
| `--entry` | | derived from CIF | Override the PDB entry ID |
| `--nf90` | | `False` | Enable NF90 mode (disables ≥ 90 % identity filter) |
| `--no-connectivity` | | connectivity on | Disable connectivity correction |
| `--tax-tsv` | | `None` | Taxonomy overrides with columns `entity_id`, `tax_id` |

*At least one of `-d` or `-m` is required.

**Outputs:** Flat gzip CSVs under `{output_dir}/`:

- `{entry_id}_seg.csv.gz` — one row per contiguous aligned range
- `{entry_id}_res.csv.gz` — one row per mapped PDB residue

---

## `db_load`

Bulk-load the per-entry CSV files produced by `segments` into DuckDB.

```bash
pdbe_sifts db_load -i INPUT_DIR -d DUCKDB
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i`, `--input-dir` | ✓ | — | Directory containing `*_seg.csv.gz` and `*_res.csv.gz` files |
| `-d`, `--duckdb` | ✓ | — | Path to the DuckDB file |

**Populates** tables `sifts_xref_segment` and `sifts_xref_residue` in the DuckDB file.

---

## `sifts2mmcif`

Inject SIFTS mappings into an mmCIF file by populating the `_pdbx_sifts_xref_db_segments` and `_pdbx_sifts_xref_db` categories.

```bash
pdbe_sifts sifts2mmcif -i CIF -o OUTPUT_DIR (-d DUCKDB | -s CSV_DIR) [options]
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i`, `--input-cif` | ✓ | — | Input CIF file |
| `-o`, `--output-dir` | ✓ | — | Output directory for SIFTS-annotated mmCIF files |
| `-d`, `--db-file` | ✗* | — | DuckDB file containing segment data |
| `-s`, `--sifts-csv-dir` | ✗* | — | Flat directory with `{entry}_seg.csv.gz` / `_res.csv.gz` files |
| `--entry` | | derived from CIF | Override the PDB entry ID |
| `-T`, `--no-track-changes` | | tracking on | Disable delta tracking against a previous run |
| `-p`, `--prev-run-dir` | | — | Directory containing a previous run's `sifts_only.mmcif` for delta comparison |

*Exactly one of `-d` or `-s` is required.

---

## `update_ncbi`

Force-update the local NCBI taxonomy database used for taxonomic scoring.

```bash
pdbe_sifts update_ncbi
```

No flags. Downloads the latest taxonomy dump from NCBI via ete4.

---

## `update_ccd_mapping`

Check the remote CCD release date and regenerate the cached three-to-one letter
mapping only when the remote source is newer.

```bash
pdbe_sifts update_ccd_mapping
```

No flags. Run `pdbe_sifts init` first so the cache path is configured.

---

## `seq2seq`

Align each polypeptide entity's canonical deposited sequence against its
coordinate sequence.

```bash
pdbe_sifts seq2seq -i CIF [-e ENTITY_ID] [-c CHAIN_ID] [--sample-fasta FASTA]
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i`, `--input-cif` | ✓ | — | Input mmCIF file |
| `-e`, `--entity-id` | | all polypeptide entities | Process one entity ID |
| `-c`, `--chain-id` | | all chains | Process one author chain ID |
| `--sample-fasta` | | `None` | Canonical sequences with headers `>{entity_id}|{chain_id1},{chain_id2}` |

---

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `PDBE_SIFTS_CONFIG` | platform config directory | Custom path to the user configuration file |
| `SIFTS_LOG_LEVEL` | `INFO` | Logging verbosity (overridden by `--log-level`) |
| `SIFTS_N_PROC` | auto | Number of internal threads per worker |
| `SIFTS_NO_CACHE_ALL` | unset | If set, disables the UniProt pickle cache |
| `SLURM_CPUS_PER_TASK` | unset | Detected automatically on SLURM clusters |
