# pdbe_sifts

Python package to run the [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) (Structure Integration with Function, Taxonomy and Sequences) pipeline locally.

Developed at [EMBL-EBI](https://www.ebi.ac.uk/) by the [PDBe](https://www.ebi.ac.uk/pdbe/) team.

---

## What is SIFTS?

SIFTS provides residue-level mappings between PDB protein structures and external databases. This package automates the core sequence alignment and mapping generation pipeline:

1. **Extract** protein sequences from PDB mmCIF files
2. **Align** them against a reference UniProt database (MMseqs2 or BLASTP)
3. **Score** and rank the alignments to select the best UniProt match per chain
4. **Generate** residue- and segment-level cross-reference mappings
5. **Export** mappings back into mmCIF files or CSV files

---

## Installation

### Recommended: Conda

```bash
conda env create -f environment.yml
conda activate pdbe_sifts
pip install -e .
```

### Manual

Install system dependencies first:
- [MMseqs2](https://github.com/soedinglab/MMseqs2) — fast sequence search
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) — optional alternative aligner
- [FASTA36] or [FASTA3]: conda install bioconda::fasta3 or (https://fasta.bioch.virginia.edu/wrpearson/fasta/)

Then install the Python package:

```bash
pip install -e .
```
or 

```bash
pip install pdbe_sifts
```

### Requirements

- Python >= 3.8
- 16 GB RAM minimum (32 GB+ recommended for large datasets)

---

## Quick Start

### 1. Initialize your config

```bash
pdbe_sifts init
```

This copies the config template to `~/.config/pdbe_sifts/config.yaml` (Linux). Edit it to set your paths.

### 2. Build a reference database

```bash
pdbe_sifts build_db \
  -i uniprot_sprot.fasta \
  -o ./my_db \
  -t taxonomy_mapping.tsv
```

The taxonomy mapping file must be tab-separated with columns `sequence_id` and `tax_id`.

### 3. Run global mappings

```bash
pdbe_sifts global_mappings \
  -i 1abc.cif \
  -od ./results \
  -db ./my_db/target_db \
  -t mmseqs \
  -threads 4
```

For a batch of entries, pass a text file listing one mmCIF path per line:

```bash
pdbe_sifts global_mappings \
  -i entries.txt \
  -od ./results \
  -db ./my_db/target_db \
  -threads 8 \
```

---

## CLI Reference

| Command | Description |
|---------|-------------|
| `pdbe_sifts init` | Copy default config to `~/.config/pdbe_sifts/config.yaml` (Linux) |
| `pdbe_sifts show` | Print the resolved configuration |
| `pdbe_sifts build_db` | Build a reference sequence database |
| `pdbe_sifts global_mappings` | Run alignment and mapping pipeline |

See [docs/cli.md](docs/cli.md) for full argument reference.

---

## Project Structure

```
src/pdbe_sifts/
├── cli.py                         # Entry point (pdbe_sifts command)
├── sifts_global_mappings.py       # Global mapping pipeline
├── sifts_segments_generation.py   # Segment generation (parallel, Batchable)
├── config/                        # OmegaConf config loading
├── base/
│   ├── batchable.py               # Parallel processing base class
│   ├── utils.py                   # Shared utilities
│   ├── log.py                     # Logging
│   └── exceptions.py              # Custom exceptions
├── mmcif/                         # mmCIF parsing (Entry, Chain, Entity, Residue)
├── global_mappings/               # Alignment tools wrappers + result parsing
│   ├── target_database.py         # Build MMseqs2/BLAST database
│   ├── mmseqs_search.py           # MMseqs2 search wrapper
│   ├── blastp.py                  # BLASTP wrapper
│   └── global_mappings_parser.py  # Parse + score alignment results
├── segments_generation/           # Alignment → CSV mapping generation
├── sifts_to_mmcif/                # Export SIFTS data to mmCIF files
├── database/
│   └── sifts_db_wrapper.py        # DuckDB wrapper (SiftsDB)
├── unp/
│   └── unp.py                     # UniProt API client + UNP class
└── data/
    └── default_config.yaml        # Default configuration template
```

---

## Documentation

- [CLI Reference](docs/cli.md)
- [Configuration](docs/configuration.md)
- [Architecture](docs/architecture.md)
- [Usage Guide](docs/usage.md)

---

## Authors

EMBL-EBI: Adam Bellaiche, Preeti Choudhary, Sreenath Sasidharan Nair, Jennifer Fleming, Sameer Velankar

## License

MIT License
