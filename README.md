# pdbe_sifts

Python package to run the [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) (Structure Integration with Function, Taxonomy and Sequences) pipeline locally.

Developed at [EMBL-EBI](https://www.ebi.ac.uk/) by the [PDBe](https://www.ebi.ac.uk/pdbe/) team.

---

## What is SIFTS?

SIFTS provides residue-level mappings between PDB protein structures and UniProt sequences. This package automates the full pipeline:

1. **Build** a reference UniProt sequence database (MMseqs2 or BLASTP)
2. **Align** PDB entity sequences against it to identify the best UniProt match per chain (≥ 90% identity) according to the SIFTS scoring function
3. **Generate** precise residue- and segment-level PDB-UNIPROT mappings via local alignment (FASTA36 `lalign36`)
4. **Store** results in a DuckDB database and per-entry CSV files
5. **Export** mappings back into annotated mmCIF files

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
```

Edit the config to set your paths (`cif_dir`, `unp_dir`, `target_db`, etc.).

### 2 — Build a reference database

```bash
pdbe_sifts build_db \
  -i uniprot_sprot.fasta \
  -o ./my_db \
  -t taxonomy_mapping.tsv   # TSV: sequence_id <tab> tax_id
```

### 3 — Run global mappings

```bash
# Single entry
pdbe_sifts global_mappings -i 1abc.cif -od ./results -db ./my_db/target_db

# Batch (one mmCIF path per line)
pdbe_sifts global_mappings -i entries.txt -od ./results -db ./my_db/target_db -threads 8
```

Produces `hits.duckdb` — a scored table of UniProt accession candidates per PDB entity.

### 4 — Generate SIFTS segments and residue mappings

```bash
# Single entry
pdbe_sifts segments -d hits.duckdb -o ./segments single --entry 1abc

# Batch (parallel, 4 workers)
pdbe_sifts segments -d hits.duckdb -o ./segments batch \
  --list entries.txt --workers 4 > segments.log 2>&1
```

Produces per-entry gzip-compressed CSV files under `{output_dir}/{entry_id}/sifts/`.

### 5 — Annotate mmCIF files with SIFTS data

```bash
pdbe_sifts sifts2mmcif -d hits.duckdb -i ./mmcif -o ./mmcif_annotated batch \
  --list entries.txt --workers 4
```

---

## CLI Reference

| Command | Description |
|---------|-------------|
| `pdbe_sifts init` | Copy default config template to `~/.config/pdbe_sifts/config.yaml` |
| `pdbe_sifts show` | Print the fully resolved configuration |
| `pdbe_sifts build_db` | Build a reference sequence database (MMseqs2 or BLASTP) from a UniProt FASTA |
| `pdbe_sifts fasta_build` | Extract entity sequences from mmCIF files and write a FASTA |
| `pdbe_sifts global_mappings` | Align PDB sequences against the reference DB; score and store hits in DuckDB |
| `pdbe_sifts segments` | Generate residue- and segment-level SIFTS mappings via local alignment (`lalign36`) |
| `pdbe_sifts sifts2mmcif` | Inject SIFTS mappings into annotated mmCIF files |

`segments` and `sifts2mmcif` accept `single` or `batch` sub-commands:

```
single --entry <pdb_id>
batch  --list <file>  --workers <N>  [--timeout <secs>]  [--no-retry]  [--failure-threshold <f>]
```

---

## Outputs

### Global mappings

| File | Format | Content |
|------|--------|---------|
| `hits.duckdb` | DuckDB | Scored UniProt accession candidates per PDB entity |
| `hits_<entry>.tsv` | TSV | Raw MMseqs2 / BLASTP alignment hits |

### Segment generation

Per entry, under `{output_dir}/{entry_id}/sifts/`:

| File | Format | Content |
|------|--------|---------|
| `sifts_segment_mapping.csv.gz` | CSV (gzip) | One row per contiguous aligned range (PDB ↔ UniProt positions, identity, conflicts, chimera flag) |
| `sifts_residue_mapping.csv.gz` | CSV (gzip) | One row per mapped PDB residue (auth seq id, UniProt position, one-letter codes, observed flag) |

When `--write-to-db` is set, results are also loaded into DuckDB tables `sifts_xref_segment` and `sifts_xref_residue`.

---

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `SIFTS_MAXTASKS_PER_CHILD` | `25` | Worker processes are recycled after this many entries to limit memory accumulation |
| `SIFTS_LALIGN_TIMEOUT` | `350` | Per-entry timeout in seconds for `lalign36` alignments |
| `SIFTS_LOG_LEVEL` | `INFO` | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`) |

---

## Project Structure

```
src/pdbe_sifts/
├── cli.py                         # CLI entry point (pdbe_sifts command)
├── sifts_global_mappings.py       # Global mapping pipeline
├── sifts_segments_generation.py   # Segment generation (SiftsAlign, Batchable subclass)
├── config/                        # OmegaConf configuration loading
├── base/
│   ├── batchable.py               # Parallel processing base class (Pool + imap_unordered)
│   ├── utils.py                   # UniProt fetch, identity/coverage helpers
│   ├── log.py                     # Logging setup
│   └── exceptions.py              # Custom exceptions
├── mmcif/                         # mmCIF parsing (Entry, Chain, Entity, Residue, ChemComp)
├── global_mappings/
│   ├── target_database.py         # Build MMseqs2 / BLAST reference database
│   ├── mmseqs_search.py           # MMseqs2 easy-search wrapper
│   ├── blastp.py                  # BLASTP wrapper
│   └── global_mappings_parser.py  # Parse TSV hits, score, store in DuckDB
├── segments_generation/
│   └── alignment/                 # lalign36 wrapper, isoform alignment, residue mapping
├── sifts_to_mmcif/                # Inject SIFTS data back into mmCIF files
├── database/
│   └── sifts_db_wrapper.py        # DuckDB schema and bulk-load helpers
├── unp/
│   └── unp.py                     # UniProt REST client, isoform handling
└── data/
    └── default_config.yaml        # Default configuration template
```

---

## Authors

EMBL-EBI PDBe team: Adam Bellaiche, Preeti Choudhary, Sreenath Sasidharan Nair, Jennifer Fleming, Sameer Velankar

## License

MIT
