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

## System Requirements

- Python >= 3.8
- Operating System: Linux, macOS, or Windows (WSL recommended)
- Minimum 4GB RAM (16GB+ recommended for large datasets)

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