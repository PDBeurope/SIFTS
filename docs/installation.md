# Installation

## System dependencies

The following external binaries must be installed and available on `PATH` before running the pipeline:

| Tool | Purpose | Recommended install |
|------|---------|---------------------|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Fast global sequence search | `conda install -c conda-forge mmseqs2` |
| [FASTA36](https://fasta.bioch.virginia.edu/wrpearson/fasta/) (`lalign36`) | Local pairwise alignment | `conda install -c bioconda fasta3` |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) (`blastp`, `makeblastdb`) | Optional alternative to MMseqs2 | `conda install -c bioconda blast` |

## Python package

=== "conda (recommended)"

    ```bash
    conda env create -f environment.yml
    conda activate pdbe_sifts
    pip install -e .
    ```

=== "pip"

    ```bash
    pip install pdbe_sifts
    ```

**Requirements:** Python ≥ 3.10 · 16 GB RAM minimum (32 GB+ recommended for large datasets)

## Post-install setup

### 1. Create your config file

```bash
pdbe_sifts init
```

This copies the built-in config template to `~/.config/pdbe_sifts/config.yaml` and downloads the NCBI taxonomy database (~70 MB, first run only).

### 2. Edit the config

Open `~/.config/pdbe_sifts/config.yaml` and set the following fields:

| Field | Description |
|-------|-------------|
| `user.base_dir` | Working directory for all pipeline outputs |
| `user.nobackup_dir` | Large-file cache directory (UniProt, CCD files) |
| `user.target_db` | Path to the pre-built reference database (after running `build_db`) |

### 3. Verify the resolved configuration

```bash
pdbe_sifts show
```

This prints the fully resolved configuration, including defaults and any overrides.
