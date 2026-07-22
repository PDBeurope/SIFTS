from importlib.metadata import PackageNotFoundError

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
        "user:\n" f"  nobackup_dir: {nobackup_dir.as_posix()}\n",
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
        "user:\n" f"  nobackup_dir: {nobackup_dir.as_posix()}\n",
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
