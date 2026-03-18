import re
from pathlib import Path

from pdbe_sifts.sifts_fasta_builder import FastaBuilder

DATA_DIR = Path(__file__).parent / "data"
EXPECTED_FASTA = DATA_DIR / "fasta" / "query.fasta"


def parse_fasta(path: Path) -> dict[str, str]:
    sequences = {}
    header = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:]
            sequences[header] = ""
        elif header is not None:
            sequences[header] += line.strip()
    return sequences


def test_build_from_cif_list(cif_list_file, tmp_path):
    fasta_path = FastaBuilder(cif_list_file, tmp_path, threads=1).build()

    assert fasta_path.exists()
    assert fasta_path.stat().st_size > 0

    produced = parse_fasta(fasta_path)
    expected = parse_fasta(EXPECTED_FASTA)
    assert produced == expected


def test_two_sequences_produced(cif_list_file, tmp_path):
    fasta_path = FastaBuilder(cif_list_file, tmp_path, threads=1).build()
    produced = parse_fasta(fasta_path)
    assert len(produced) == 2


def test_fasta_header_format(cif_list_file, tmp_path):
    fasta_path = FastaBuilder(cif_list_file, tmp_path, threads=1).build()
    produced = parse_fasta(fasta_path)
    pattern = re.compile(r"^pdb\|[a-z0-9]+-\d+\|OX=\d+$")
    for header in produced:
        assert pattern.match(
            header
        ), f"Header does not match expected format: {header}"


def test_fasta_passthrough(tmp_path):
    fasta_path = FastaBuilder(EXPECTED_FASTA, tmp_path).build()
    assert fasta_path == EXPECTED_FASTA
