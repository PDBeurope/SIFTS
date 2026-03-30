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

### System dependencies

The following binaries must be installed and available on `PATH`:

| Tool | Purpose | Install |
|------|---------|---------|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Fast global sequence search | `conda install -c conda-forge mmseqs2` |
| [FASTA36](https://fasta.bioch.virginia.edu/wrpearson/fasta/) (`lalign36`) | Local pairwise alignment | `conda install -c bioconda fasta3` |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | Optional alternative to MMseqs2 | `conda install -c bioconda blast` |

### Python package

```bash
# Recommended: conda environment
conda env create -f environment.yml
conda activate pdbe_sifts
pip install -e .

# Or directly
pip install pdbe_sifts
```

**Requirements:** Python ≥ 3.10 · 16 GB RAM minimum (32 GB+ recommended for large datasets)

---

## Quick Start

### 1 — Initialise your config

```bash
pdbe_sifts init
# → creates ~/.config/pdbe_sifts/config.yaml
# → downloads the NCBI taxonomy database (~70 MB, first run only)
```

Edit the config to set your paths (`base_dir`, `nobackup_dir`, `target_db` (after building it), etc.). You can also setup different alignment parameters.

### 2 — Build a reference database

```bash
pdbe_sifts build_db \
  -i uniprot_sprot.fasta \
  -o ./my_db \
  -t taxonomy_mapping.tsv   # TSV: sequence_id <tab> tax_id
```

### 3 — Run structure to sequence matching

```bash
# Single CIF entry
pdbe_sifts sequence_match -i 1abc.cif -o ./results -d ./my_db/target_db

# Batch (one mmCIF path per line)
pdbe_sifts sequence_match -i entries.txt -o ./results -d ./my_db/target_db --threads 8
```

At this step you can also provide a .csv file to faster the scoring function. This CSV file must contains per row: row_num, uniprot_accession, dataset (Swiss-Prot or TrEMBL), pdb cross-references, annotation score.

Produces `hits.duckdb` and `hits.tsv` — a scored and raw table of sequence candidates per structure entity.

### 4 — Generate SIFTS segments and residue mappings

```bash
# With DuckDB hits (from structure to sequence matching step)
pdbe_sifts segments -i 1abc.cif.gz -o ./segments -d hits.duckdb

# Manual structure-sequence mapping (chain:accession)
pdbe_sifts segments -i 1abc.cif.gz -o ./segments -m "A:P00963,B:P00963"

# Custom FASTA mapping (headers: >{structure_id}|{auth_asym_id}|{sequence_id})
pdbe_sifts segments -i 1abc.cif.gz -o ./segments -m custom_seqs.fasta
```

Produces per-entry gzip-compressed CSV files under `{output_dir}/`.

### 5 — Load segment data into DuckDB

```bash
pdbe_sifts db_load -i ./segments/ -d hits.duckdb
```

Bulk-loads the segment and residue CSVs produced in step 4 into the `sifts_xref_segment` and `sifts_xref_residue` tables of the DuckDB file.

### 6 — Annotate mmCIF files with residue level mappings and SIFTS data

```bash
# Reading from DuckDB (after step 5)
pdbe_sifts sifts2mmcif \
  -i 1abc.cif.gz \
  -o ./sifts_mmcif \
  -d hits.duckdb

# Or reading segment CSVs directly (skip step 5)
pdbe_sifts sifts2mmcif \
  -i 1abc.cif.gz \
  -o ./sifts_mmcif \
  -s ./segments/
```

---

## CLI Reference

| Command | Description |
|---------|-------------|
| `pdbe_sifts init` | Copy default config to `~/.config/pdbe_sifts/config.yaml` and init NCBI taxonomy DB |
| `pdbe_sifts show` | Print the fully resolved configuration |
| `pdbe_sifts update_ncbi` | Force-update the local NCBI taxonomy database (ete4) |
| `pdbe_sifts build_db` | Build a reference sequence database (MMseqs2 or BLASTP) from a FASTA file |
| `pdbe_sifts fasta_build` | Extract entity sequences from mmCIF files and write a FASTA |
| `pdbe_sifts sequence_match` | Align structure sequences against the reference DB; score and store hits in DuckDB |
| `pdbe_sifts segments` | Generate SIFTS mappings for a **single** mmCIF entry |
| `pdbe_sifts db_load` | Bulk-load segment/residue CSVs from segments generation into DuckDB |
| `pdbe_sifts sifts2mmcif` | Inject SIFTS mappings into an annotated mmCIF file |

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
# → writes hits.duckdb and hits_<entry>.tsv to out_dir
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

# Mode 2: provide a manual mapping (accessions or custom FASTA)
sa = SiftsAlign(
    cif_file="1abc.cif",
    out_dir="./segments/",
    unp_mode="A:P00963,B:P00963",   # or path to a FASTA file
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
├── config/                        # OmegaConf configuration loading
├── base/
│   ├── paths.py                   # Single load_config() + all configuration getters
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

MIT
