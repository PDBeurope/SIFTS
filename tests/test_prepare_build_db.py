import gzip

import pytest
from pdbe_sifts.prepare_build_db import (
    parse_uniprot_tax_header,
    prepare_build_db,
    write_uniprot_tax_mapping,
)


def test_parse_uniprot_tax_header_supports_sp_and_tr():
    assert parse_uniprot_tax_header(
        "sp|P29373|RABP2_HUMAN Protein OS=Homo sapiens OX=9606 PE=1"
    ) == ("P29373", "9606")
    assert parse_uniprot_tax_header(
        "tr|A0A000|A0A000_HUMAN Protein OS=Homo sapiens OX=9606 PE=1"
    ) == ("A0A000", "9606")


def test_prepare_build_db_with_input_fasta_creates_tax_mapping(tmp_path):
    input_fasta = tmp_path / "input.fasta.gz"
    output_fasta = tmp_path / "prepared.fasta.gz"
    output_tax_mapping = tmp_path / "tax_mapping.tsv"

    with gzip.open(input_fasta, "wt", encoding="utf-8") as handle:
        handle.write(
            ">sp|P83346|3SO3_BUNCA Bucain OS=Bungarus candidus OX=92438 PE=1 SV=1\n"
            "RKCLIKYS\n"
            ">tr|A0A000|A0A000_HUMAN Protein OS=Homo sapiens OX=9606 PE=1 SV=1\n"
            "MPEPTIDE\n"
            ">xx|SKIPME|Unsupported OS=Example OX=1\n"
            "M\n"
        )

    fasta_path, tax_mapping_path = prepare_build_db(
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        output_tax_mapping=output_tax_mapping,
    )

    assert fasta_path == output_fasta
    assert tax_mapping_path == output_tax_mapping
    assert output_fasta.read_bytes() == input_fasta.read_bytes()
    assert output_tax_mapping.read_text(encoding="utf-8").splitlines() == [
        "P83346\t92438",
        "A0A000\t9606",
    ]


def test_prepare_build_db_compresses_plain_input_when_output_is_gzip(tmp_path):
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "prepared.fasta.gz"
    output_tax_mapping = tmp_path / "tax_mapping.tsv"
    input_fasta.write_text(
        ">sp|P29373|RABP2_HUMAN Protein OS=Homo sapiens OX=9606 PE=1\n"
        "MPEPTIDE\n",
        encoding="utf-8",
    )

    prepare_build_db(
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        output_tax_mapping=output_tax_mapping,
    )

    with gzip.open(output_fasta, "rt", encoding="utf-8") as handle:
        assert handle.readline().startswith(">sp|P29373|")
    assert output_tax_mapping.read_text(encoding="utf-8") == "P29373\t9606\n"


@pytest.mark.parametrize("compressed", [False, True])
def test_write_uniprot_tax_mapping_supports_plain_and_gzip_fasta(
    tmp_path, compressed
):
    suffix = ".fasta.gz" if compressed else ".fasta"
    input_fasta = tmp_path / f"input{suffix}"
    output_tax_mapping = tmp_path / "nested" / "tax_mapping.tsv"
    content = (
        ">sp|P83346|3SO3_BUNCA Bucain OS=Bungarus candidus OX=92438 PE=1\n"
        "RKCLIKYS\n"
        ">xx|SKIPME|Unsupported OS=Example OX=1\n"
        "M\n"
        ">tr|A0A000|A0A000_HUMAN Protein OS=Homo sapiens OX=9606 PE=1\n"
        "MPEPTIDE\n"
    )
    if compressed:
        with gzip.open(input_fasta, "wt", encoding="utf-8") as handle:
            handle.write(content)
    else:
        input_fasta.write_text(content, encoding="utf-8")

    result = write_uniprot_tax_mapping(input_fasta, output_tax_mapping)

    assert result == output_tax_mapping
    assert output_tax_mapping.read_text(encoding="utf-8") == (
        "P83346\t92438\nA0A000\t9606\n"
    )


def test_write_uniprot_tax_mapping_rejects_fasta_without_taxonomy_ids(
    tmp_path,
):
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">custom_identifier\nMPEPTIDE\n", encoding="utf-8")

    with pytest.raises(
        ValueError, match="No UniProtKB taxonomy IDs were parsed"
    ):
        write_uniprot_tax_mapping(input_fasta, tmp_path / "tax_mapping.tsv")
