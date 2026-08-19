# Installation

## System dependencies

The following external binaries must be installed and available on `PATH` before running the pipeline:

| Tool | Purpose | Recommended install |
|------|---------|---------------------|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Fast global sequence search | `conda install -c conda-forge mmseqs2` |
| [FASTA36](https://fasta.bioch.virginia.edu/wrpearson/fasta/) (`lalign36`) | Local pairwise alignment | `conda install -c bioconda fasta3` |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) (`blastp`, `makeblastdb`) | Optional alternative to MMseqs2 | `conda install -c bioconda blast` |

## Python package

=== "Local checkout with conda (recommended)"

    ```bash
    git clone https://github.com/PDBeurope/SIFTS
    cd SIFTS
    conda env create -f environment.yml
    conda activate pdbe_sifts

    # Replace the published package installed by environment.yml with
    # this checkout in editable mode.
    pip install -e .
    ```

=== "Published package"

    ```bash
    pip install pdbe-sifts
    ```

=== "Local checkout with uv"

    ```bash
    git clone https://github.com/PDBeurope/SIFTS
    cd SIFTS
    uv sync
    source .venv/bin/activate  # macOS/Linux
    ```

    `uv sync` creates `.venv` and `uv.lock`, installs the dependencies,
    and installs the checkout in editable mode.

The external binaries listed above are not installed by `pip` or `uv`. BLAST+
is optional unless `--tool blastp` is used.

**Requirements:** Python ≥ 3.10 · 16 GB RAM minimum (32 GB+ recommended for large datasets)

## Manual installation of external binaries

BLAST+:

```bash
brew install blast          # macOS
sudo apt install ncbi-blast+ # Debian/Ubuntu
```

MMseqs2:

```bash
brew install mmseqs2        # macOS or Linux with Homebrew
```

FASTA36 from source (Linux):

```bash
git clone https://github.com/wrpearson/fasta36.git
cd fasta36/src
make -f ../make/Makefile.linux_sse2 all
export PATH="$(pwd)/../bin:$PATH"
```

## Verify the installation

```bash
pdbe_sifts --version
command -v mmseqs
command -v lalign36

# Required only when using --tool blastp
command -v makeblastdb
```

## Post-install setup

### 1. Create your config file

```bash
pdbe_sifts init
```

This copies the built-in config template to
`~/.config/pdbe_sifts/config.yaml`, initialises the NCBI taxonomy database,
builds the UniProt–PDB DuckDB index, and generates the CCD mapping cache.

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
