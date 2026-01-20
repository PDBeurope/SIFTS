# pdbe_sifts

Python package to run SIFTS from EMBL-EBI (PDBe) locally.

## Installation

### System Requirements

This package requires several bioinformatics tools that cannot be installed via pip alone. We strongly recommend using conda/mamba for installation.

### Installation with Conda (Recommended)

1. Create the conda environment from the YAML file:
```bash
conda env create -f environment.yml
conda activate pdbe_sifts
```

2. Install the package in development mode:
```bash
pip install -e .
```

### Manual Installation (Alternative)

If you prefer not to use conda, you must manually install the following system dependencies:

#### Required System Dependencies:
- **MMseqs2**: Sequence similarity search tool
- **BLAST**: Basic Local Alignment Search Tool
- **FASTA3**: Sequence comparison tool suite

#### Installation on Linux (Ubuntu/Debian):
```bash
# MMseqs2
# Download from: https://github.com/soedinglab/MMseqs2/releases
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
tar xvfz mmseqs-linux-avx2.tar.gz
export PATH=$(pwd)/mmseqs/bin/:$PATH

# BLAST
sudo apt-get install ncbi-blast+

# FASTA3
sudo apt-get install fasta3
```

#### Installation on macOS:
```bash
# With Homebrew
brew install mmseqs2
brew install blast
brew install brewsci/bio/fasta3
```

#### Then install Python dependencies:
```bash
pip install -r requirements.txt
pip install -e .
```

## Python Dependencies

The following dependencies will be automatically installed via pip:

- tqdm
- biopython
- funcy
- lxml
- requests
- pandas
- pymmseqs
- coloredlogs
- pyyaml
- scikit-learn
- ete4
- omegaconf
- gemmi

## Verifying Installation

To verify that all tools are correctly installed:
```bash
python -c "import pdbe_sifts; print('Package imported successfully')"
mmseqs --version
blastp -version
fasta36 -version
```

## System Requirements

- Python >= 3.8
- Operating System: Linux, macOS, or Windows (WSL recommended)
- Minimum 4GB RAM (8GB+ recommended for large datasets)

## Usage

[Add basic usage examples here]
```python
import pdbe_sifts

# Example usage
```

## Support

For questions or issues, please contact: adamb@ebi.ac.uk

## License

MIT License

## Authors

EMBL-EBI: Adam Bellaiche, Preeti Choudhary, Sreenath Sasidharan Nair, Jennifer Fleming, Sameer Velankar

## Citation

[Add citation information if applicable]