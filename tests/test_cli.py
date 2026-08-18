import sys
from importlib.metadata import PackageNotFoundError
from unittest.mock import MagicMock

import pytest
from pdbe_sifts import cli


def test_get_pdbe_sifts_version_from_metadata(monkeypatch):
    monkeypatch.setattr(cli, "metadata_version", lambda name: "9.9.9")

    assert cli._get_pdbe_sifts_version() == "9.9.9"


def test_get_pdbe_sifts_version_falls_back_to_package_version(monkeypatch):
    def raise_not_found(name):
        raise PackageNotFoundError(name)

    monkeypatch.setattr(cli, "metadata_version", raise_not_found)
    monkeypatch.setattr(cli, "__version__", "1.2.3")

    assert cli._get_pdbe_sifts_version() == "1.2.3"


def test_setup_cache_creates_configured_directories(tmp_path):
    config_path = tmp_path / "config.yaml"
    nobackup_dir = tmp_path / "nobackup"
    config_path.write_text(
        f"user:\n  nobackup_dir: {nobackup_dir.as_posix()}\n",
        encoding="utf-8",
    )

    cache_paths = cli._setup_cache(config_path)

    assert cache_paths["base"] == nobackup_dir / "sifts_data_cache"
    assert cache_paths["uniprot"].is_dir()
    assert cache_paths["ccd"].is_dir()
    assert cache_paths["three_to_one"].parent.is_dir()


def test_setup_cache_no_create_only_resolves_paths(tmp_path):
    config_path = tmp_path / "config.yaml"
    nobackup_dir = tmp_path / "nobackup"
    config_path.write_text(
        f"user:\n  nobackup_dir: {nobackup_dir.as_posix()}\n",
        encoding="utf-8",
    )

    cache_paths = cli._setup_cache(config_path, create=False)

    assert cache_paths["base"] == nobackup_dir / "sifts_data_cache"
    assert not cache_paths["base"].exists()


def test_setup_cache_rejects_unconfigured_cache(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text("user:\n  nobackup_dir: null\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Cache configuration is incomplete"):
        cli._setup_cache(config_path)


def test_mapping_fasta_path_rejects_legacy_literal():
    with pytest.raises(
        cli.argparse.ArgumentTypeError,
        match="mapping FASTA file does not exist",
    ):
        cli._mapping_fasta_path("A:P29373")


def test_segments_cli_rejects_legacy_mapping(monkeypatch, capsys):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pdbe_sifts",
            "segments",
            "-i",
            "entry.cif",
            "-o",
            "segments",
            "-m",
            "A:P29373",
        ],
    )

    with pytest.raises(SystemExit, match="2"):
        cli.main()

    assert (
        "mapping FASTA file does not exist: A:P29373" in capsys.readouterr().err
    )


def test_segments_cli_passes_mapping_fasta(tmp_path, monkeypatch):
    cif_path = tmp_path / "entry.cif"
    cif_path.write_text("data_1cbs\n_entry.id 1cbs\n", encoding="utf-8")
    fasta_path = tmp_path / "mapping.fasta"
    fasta_path.write_text(">1cbs|A|MYREF\nACDE\n", encoding="utf-8")
    output_dir = tmp_path / "segments"
    align = MagicMock()
    align.conn = None
    constructor = MagicMock(return_value=align)
    monkeypatch.setattr(
        "pdbe_sifts.sifts_segments_generation.SiftsAlign", constructor
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pdbe_sifts",
            "segments",
            "-i",
            str(cif_path),
            "-o",
            str(output_dir),
            "--entry",
            "1cbs",
            "-m",
            str(fasta_path),
        ],
    )

    cli.main()

    assert constructor.call_args.kwargs["mapping_fasta"] == fasta_path
    align.process_entry.assert_called_once_with("1cbs")


def test_create_tax_file_cli_writes_mapping_and_reports_path(
    tmp_path, monkeypatch, capsys
):
    input_fasta = tmp_path / "input.fasta"
    output_tax_mapping = tmp_path / "taxonomy" / "mapping.tsv"
    input_fasta.write_text(
        ">sp|P29373|RABP2_HUMAN Protein OS=Homo sapiens OX=9606 PE=1\n"
        "MPEPTIDE\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pdbe_sifts",
            "create_tax_file",
            "--input-fasta",
            str(input_fasta),
            "--output-tax-mapping",
            str(output_tax_mapping),
        ],
    )

    cli.main()

    assert output_tax_mapping.read_text(encoding="utf-8") == "P29373\t9606\n"
    assert capsys.readouterr().out == (
        f"Taxonomy mapping written to: {output_tax_mapping}\n"
    )
