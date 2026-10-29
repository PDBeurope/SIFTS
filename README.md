# pdbe_sifts

Python package to run the [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) (Structure Integration with Function, Taxonomy and Sequences) pipeline locally.

Developed at [EMBL-EBI](https://www.ebi.ac.uk/) by the [PDBe](https://www.ebi.ac.uk/pdbe/) team.

---

## What is SIFTS?

SIFTS provides residue-level mappings between structures and sequences. This package automates the full pipeline:

1. **Build** a reference sequence database (MMseqs2 or BLASTP)
2. **Align** structure sequences against it to identify the best match per chain (≥ 90% identity) according to the SIFTS scoring function
3. **Generate** precise residue- and segment-level structure-sequence mappings via local alignment (FASTA36 `lalign36`)
4. **Store** results in a DuckDB database and per-entry CSV files
5. **Export** mappings back into annotated mmCIF files

The whole pipeline can work on non-UniProt or non-PDB entries. However, it will use only the adjusted score to rank the hits.

---

## Installation

**Requirements:** Python ≥ 3.10 · 16 GB RAM minimum (32 GB+ recommended for large datasets)

Besides the Python package, the pipeline calls these external tools, which
must be on `PATH`:

| Tool | Purpose |
|------|---------|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) (`mmseqs`) | Fast global sequence search |
| [FASTA36](https://fasta.bioch.virginia.edu/wrpearson/fasta/) (`lalign36`) | Local pairwise alignment |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) (`blastp`, `makeblastdb`) | Optional alternative to MMseqs2, used with `--tool blastp` |

The simplest way to get everything is a
[micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
environment with the tools from conda-forge and bioconda:

```bash
micromamba create -n sifts -c conda-forge -c bioconda python=3.13 pip mmseqs2 fasta3 blast
micromamba activate sifts
pip install pdbe-sifts
```

`conda` and `mamba` accept the same arguments if you already use one of them.
If you install the tools another way (see the
[installation guide](docs/installation.md)), `pip install pdbe-sifts` is all
that remains. Windows is not supported; use the Windows Subsystem for Linux.

## Development setup

The development environment is pinned in `conda-lock.yml` for linux-64 and
osx-arm64. It contains the runtime dependencies, the external tools, and the
test, lint and docs tooling.

```bash
git clone https://github.com/PDBeurope/SIFTS
cd SIFTS

micromamba create -n sifts-dev -f conda-lock.yml --category "main dev docs"
micromamba activate sifts-dev

pip install --no-deps --no-build-isolation -e .
pre-commit install
```

Run the tests and preview the docs:

```bash
pytest
mkdocs serve
```

To change a dependency, edit `pyproject.toml` (the Python packages are
listed under `[project]`, the conda-only tools and Python version under
`[tool.conda-lock.dependencies]`), then refresh the lock file:

```bash
conda-lock lock --micromamba -f pyproject.toml -e dev -e docs \
  --pypi_to_conda_lookup_file https://raw.githubusercontent.com/conda-forge/conda-forge-bot-data/main/mappings/pypi/grayskull_pypi_mapping.json
```

The `conda-lock` pre-commit hook runs the same command, so a commit that
changes `pyproject.toml` without updating `conda-lock.yml` fails.

## Quick Start

### 1.1 — Initialise your required files to run PDBe-SIFTS

You first need a fasta file containing sequences and a taxonomic mapping file. For that you can use the command:

```bash
pdbe_sifts prepare_build_db --output-fasta /my/path/to/myfile.fasta.gz --output-tax-mapping /my/path/to/taxonomy_mapping.tsv
```

You can also specify an existing FASTA file using `--input-fasta`. Without it,
the command downloads the UniProtKB/Swiss-Prot FASTA file by default.

This command assumes that the FASTA headers follow UniProtKB conventions.

If the FASTA is already in its final location and only the taxonomy mapping is
needed, create it without copying or downloading the FASTA:

```bash
pdbe_sifts create_tax_file \
  --input-fasta /my/path/to/myfile.fasta.gz \
  --output-tax-mapping /my/path/to/taxonomy_mapping.tsv
```

Both plain and gzip-compressed FASTA files are accepted. Headers must use the
UniProtKB `sp|...` or `tr|...` format and contain `OX=<taxid>`. The output has
no header and contains one `ACCESSION<TAB>TAXID` row per compatible header.

### 1.2 — Initialise your config

```bash
pdbe_sifts init
# → creates ~/.sifts/config.yaml
# → downloads the NCBI taxonomy database (~70 MB, first run only)
# → builds the UniProt–PDB cross-reference index and CCD mapping cache
```

By default, `init` sets `base_dir` to the current working directory, `nobackup_dir` to `<base_dir>/cache` and `target_db` to `<base_dir>/database/sp.db`. Override them with `--base-dir`, `--nobackup-dir` and `--target-db`, or edit the config afterwards. You can also setup different alignment parameters.

`base_dir`: path to where you want the results to be saved

`nobackup_dir`:path to where you want the tmp files to be saved

`target_db`: path to the database created using build_db

### 2 — Build a reference database

If you already have a UniProtKB FASTA file, the complete path from that FASTA
to a searchable database is:

```bash
pdbe_sifts create_tax_file \
  --input-fasta /my/path/to/myfile.fasta.gz \
  --output-tax-mapping /my/path/to/taxonomy_mapping.tsv

pdbe_sifts build_db \
  -i /my/path/to/myfile.fasta.gz \
  -t /my/path/to/taxonomy_mapping.tsv \
  --tool mmseqs \
  --threads 8
```

`create_tax_file` requires UniProtKB headers in the format described above. For
a custom FASTA header format, create the headerless taxonomy TSV yourself, with
one `sequence_id<TAB>tax_id` row per sequence, and pass it directly to
`build_db` with `-t`.

`build_db` accepts a gzip-compressed FASTA when using MMseqs2. It creates
several database files using `user.target_db` (by default
`<base_dir>/database/sp.db`) as their path prefix; pass `-o` to override it. MMseqs2 is
the default backend; use `--tool blastp` to build a BLAST database instead.

### 3 — Run structure to sequence matching

```bash
# Single CIF entry
pdbe_sifts sequence_match -i 1abc.cif -o ./results -d ./my/db/folder/prefix

# Batch (one mmCIF path per line)
pdbe_sifts sequence_match -i entries.txt -o ./results -d ./my/db/folder/prefix --threads 8
```

At this step you can also provide a .csv file to faster the scoring function. This CSV file must contains per row: row_num, uniprot_accession, dataset (Swiss-Prot or TrEMBL), pdb cross-references, annotation score.

This command produces `hits.duckdb` and `hits_<entry>.tsv` — scored and raw
sequence candidates per structure entity. The files are written under
`{output_dir}/{tool}_{input_name}/`; for example,
`./results/mmseqs_entries/hits.duckdb`.

### 4 — Generate SIFTS segments and residue mappings

```bash
# With DuckDB hits (from structure to sequence matching step)
pdbe_sifts segments \
  -i 1abc.cif.gz \
  -o ./segments \
  -d ./results/mmseqs_1abc/hits.duckdb

# Custom FASTA mapping (headers: >{structure_id}|{auth_asym_id}|{sequence_id})
pdbe_sifts segments -i 1abc.cif.gz -o ./segments -m custom_seqs.fasta
```

Produces per-entry gzip-compressed CSV files under `{output_dir}/`.

### 5 — Load segment data into DuckDB

```bash
pdbe_sifts db_load \
  -i ./segments/ \
  -d ./results/mmseqs_1abc/hits.duckdb
```

Bulk-loads the segment and residue CSVs produced in step 4 into the `sifts_xref_segment` and `sifts_xref_residue` tables of the DuckDB file.

### 6 — Annotate mmCIF files with residue level mappings and SIFTS data

If you run the `db_load` command, you can use the command below because the `-d` argument requires two tables created by the db_load
step.

```bash
# Reading from DuckDB (after step 5)
pdbe_sifts sifts2mmcif \
  -i 1abc.cif.gz \
  -o ./sifts_mmcif \
  -d ./results/mmseqs_1abc/hits.duckdb
```

If not please use:
```bash
# Or reading segment CSVs directly (skip step 5)
pdbe_sifts sifts2mmcif \
  -i 1abc.cif.gz \
  -o ./sifts_mmcif \
  -s ./segments/
```

This integrates the computed SIFTS annotations into the source PDBx/mmCIF file by populating the `_pdbx_sifts_unp_segments`. Additionally, the `atom_site` category is extended to include the sequence accession and residue number from the best available mapping, enabling each atomic coordinate to be associated directly with the corresponding mapped residue.

---

## CLI Reference

| Command | Description |
|---------|-------------|
| `pdbe_sifts init` | Create the config and initialise the taxonomy, UniProt–PDB, and CCD caches |
| `pdbe_sifts show` | Print the fully resolved configuration |
| `pdbe_sifts update_ncbi` | Force-update the local NCBI taxonomy database (ete4) |
| `pdbe_sifts prepare_build_db` | Prepare FASTA and taxonomy mapping inputs for target database creation |
| `pdbe_sifts create_tax_file` | Extract an accession-to-taxid TSV from a UniProtKB FASTA file |
| `pdbe_sifts build_db` | Build a reference sequence database (MMseqs2 or BLASTP) from a FASTA file |
| `pdbe_sifts fasta_build` | Extract entity sequences from mmCIF files and write a FASTA |
| `pdbe_sifts sequence_match` | Align structure sequences against the reference DB; score and store hits in DuckDB |
| `pdbe_sifts segments` | Generate SIFTS mappings for a **single** mmCIF entry |
| `pdbe_sifts db_load` | Bulk-load segment/residue CSVs from segments generation into DuckDB |
| `pdbe_sifts sifts2mmcif` | Inject SIFTS mappings into an annotated mmCIF file |
| `pdbe_sifts update_ccd_mapping` | Check whether the remote CCD file is newer than the cached three-to-one letter mapping CSV and regenerate it if so |
| `pdbe_sifts seq2seq` | Align canonical deposited sequence vs coordinate sequence |

---

## Useful Classes

The pipeline classes can be used directly in Python scripts without going through the CLI.

### `TargetDb` — Build a reference sequence database

```python
from pdbe_sifts.sequence_match.target_database import TargetDb

TargetDb(
    input_path="uniprot_sprot.fasta",
    output_path="./my_db/target_db",
    tax_mapping_file="taxonomy.tsv",
    tool="mmseqs",   # or "blastp"
    threads=8,
).run()
```

### `FastaBuilder` — Extract sequences from mmCIF files

```python
from pdbe_sifts.sifts_fasta_builder import FastaBuilder

fasta_path = FastaBuilder(
    input_path="1abc.cif",   # or .cif.gz, or a .txt file listing CIF paths
    out_dir="./fasta/",
    threads=4,
).build()
```

### `SiftsSequenceMatch` — Run the alignment and scoring pipeline

```python
from pdbe_sifts.sifts_sequence_match import SiftsSequenceMatch

SiftsSequenceMatch(
    input_file="1abc.cif",   # or .fasta, or a .txt list of CIF paths
    out_dir="./results/",
    db_file="./my_db/target_db",
    tool="mmseqs",           # or "blastp"
    threads=8,
).process()
# → writes hits.duckdb and hits_<entry>.tsv under out_dir/mmseqs_1abc/
```

### `SiftsAlign` — Generate per-entry segment and residue mappings

```python
from pdbe_sifts.sifts_segments_generation import SiftsAlign

# Mode 1: use scored hits from sequence_match
sa = SiftsAlign(
    cif_file="1abc.cif",
    out_dir="./segments/",
    db_conn_str="hits.duckdb",
)

# Mode 2: provide a custom FASTA mapping
sa = SiftsAlign(
    cif_file="1abc.cif",
    out_dir="./segments/",
    mapping_fasta="custom_seqs.fasta",
)

sa.process_entry("1abc")
if sa.conn:
    sa.conn.close()
# → writes {out_dir}/1abc_seg.csv.gz
#           {out_dir}/1abc_res.csv.gz
```

### `SiftsDB` — Bulk-load segment CSVs into DuckDB

```python
import duckdb
from pdbe_sifts.database.sifts_db_wrapper import SiftsDB

conn = duckdb.connect("hits.duckdb")
SiftsDB(conn).bulk_load_from_entries("./segments/")
conn.close()
```

---

## Outputs

### Global mappings

| File | Format | Content |
|------|--------|---------|
| `hits.duckdb` | DuckDB | Scored sequence accession candidates per structure entity |
| `hits_<entry>.tsv` | TSV | Raw MMseqs2 / BLASTP alignment hits |

### Segment generation

Per entry, under `{output_dir}`:

| File | Format | Content |
|------|--------|---------|
| `{entry}_seg.csv.gz` | CSV (gzip) | One row per contiguous aligned range (structure ↔ sequence positions, identity, conflicts, chimera flag) |
| `{entry}_res.csv.gz` | CSV (gzip) | One row per mapped structure residue (auth seq id, sequence position, one-letter codes, observed flag) |
| `{entry}_nf90_seg.csv.gz` | CSV (gzip) | NF90 variant of the segment file (written when applicable) |

After running `db_load`, results are available in DuckDB tables `sifts_xref_segment` and `sifts_xref_residue`.

---

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `PDBE_SIFTS_CONFIG` | `~/.sifts/config.yaml` | Custom path to the user configuration file |
| `SIFTS_LOG_LEVEL` | `INFO` | Logging verbosity: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL` |
| `SIFTS_N_PROC` | auto | Number of internal threads per worker (lalign36 jobs). Override manually to cap CPU use. |
| `SIFTS_NO_CACHE_ALL` | unset | If set (any value), disables the UniProt pickle cache and always fetches from the REST API. |
| `SLURM_CPUS_PER_TASK` | unset | Detected automatically on SLURM clusters. Used by `get_allocated_cpus()` to set the thread count when running under a SLURM job allocation. |

---

## Project Structure

```
src/pdbe_sifts/
├── cli.py                         # CLI entry point (pdbe_sifts command)
├── sifts_sequence_match.py       # Global mapping pipeline (SiftsSequenceMatch)
├── sifts_segments_generation.py   # Single-entry segment generation (SiftsAlign)
├── sifts_fasta_builder.py         # Extract sequences from mmCIF → FASTA (FastaBuilder)
├── sifts_database_loader.py       # Standalone bulk-loader script (wraps SiftsDB)
├── config/                        # OmegaConf configuration loading — defines load_config()
├── base/
│   ├── paths.py                   # All configuration getters (imports load_config() from config/)
│   ├── utils.py                   # UniProt fetch, CPU helpers, SiftsAction
│   ├── log.py                     # Logging setup (StreamHandler, coloredlogs)
│   └── exceptions.py              # All custom exceptions (centralised)
├── database/
│   └── sifts_db_wrapper.py        # SiftsDB: DuckDB schema + bulk loader
├── mmcif/                         # mmCIF parsing (Entry, Chain, Entity, Residue, ChemComp)
├── sequence_match/
│   ├── target_database.py         # Build MMseqs2 / BLAST reference database (TargetDb)
│   ├── mmseqs_search.py           # MMseqs2 easy-search wrapper
│   ├── blastp.py                  # BLASTP wrapper
│   └── sequence_match_parser.py  # Parse TSV hits, score, store in DuckDB
├── segments_generation/
│   └── alignment/                 # lalign36 wrapper, isoform alignment, residue mapping
├── sifts_to_mmcif/                # Inject SIFTS data back into mmCIF files
├── unp/
│   └── unp.py                     # UniProt REST client, pickle cache, isoform handling
└── data/
    └── default_config.yaml        # Default configuration template (all tuneable params)
```

---

## Authors

EMBL-EBI PDBe team: Adam Bellaiche, Preeti Choudhary, Sreenath Sasidharan Nair, Jennifer Fleming, Sameer Velankar

## License

Apache-2.0
