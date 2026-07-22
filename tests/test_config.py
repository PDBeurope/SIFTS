from importlib import reload
from pathlib import Path

import pytest
from pdbe_sifts.config import init_config, load_config, set_unp_pdb_xrefs_path


def write_config(path: Path, base_dir: Path) -> None:
    path.write_text(
        "user:\n"
        f"  base_dir: {base_dir.as_posix()}\n"
        f"  nobackup_dir: {(base_dir / 'cache').as_posix()}\n",
        encoding="utf-8",
    )


def test_load_config_uses_pdbe_sifts_config_env_var(tmp_path, monkeypatch):
    config_path = tmp_path / "env_config.yaml"
    base_dir = tmp_path / "env_base"
    write_config(config_path, base_dir)
    monkeypatch.setenv("PDBE_SIFTS_CONFIG", str(config_path))

    cfg = load_config()

    assert cfg.user.base_dir == str(base_dir)
    assert cfg.user.nobackup_dir == str(base_dir / "cache")


def test_init_config_uses_pdbe_sifts_config_env_var(tmp_path, monkeypatch):
    config_path = tmp_path / "pdbe_sifts_config.yaml"
    xref_db = tmp_path / "uniprot_pdb.duckdb"
    monkeypatch.setenv("PDBE_SIFTS_CONFIG", str(config_path))

    written_path = init_config()
    set_unp_pdb_xrefs_path(xref_db)
    cfg = load_config()

    assert written_path == config_path
    assert config_path.exists()
    assert cfg.user.unp_pdb_xrefs == str(xref_db)


def test_load_config_explicit_path_takes_priority_over_env_var(
    tmp_path, monkeypatch
):
    env_config = tmp_path / "env_config.yaml"
    explicit_config = tmp_path / "explicit_config.yaml"
    write_config(env_config, tmp_path / "env_base")
    write_config(explicit_config, tmp_path / "explicit_base")
    monkeypatch.setenv("PDBE_SIFTS_CONFIG", str(env_config))

    cfg = load_config(explicit_config)

    assert cfg.user.base_dir == str(tmp_path / "explicit_base")


def test_load_config_rejects_missing_pdbe_sifts_config(monkeypatch):
    missing_config = Path("/does/not/exist/pdbe_sifts_config.yaml")
    monkeypatch.setenv("PDBE_SIFTS_CONFIG", str(missing_config))

    with pytest.raises(
        FileNotFoundError,
        match="PDBE_SIFTS_CONFIG points to missing file",
    ):
        load_config()


def test_paths_module_uses_pdbe_sifts_config_env_var(tmp_path, monkeypatch):
    config_path = tmp_path / "env_config.yaml"
    base_dir = tmp_path / "env_base"
    write_config(config_path, base_dir)
    monkeypatch.setenv("PDBE_SIFTS_CONFIG", str(config_path))

    import pdbe_sifts.base.paths as paths

    reload(paths)

    assert paths.get_conf_user_base_dir() == str(base_dir)

    monkeypatch.delenv("PDBE_SIFTS_CONFIG")
    reload(paths)
