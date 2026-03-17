from __future__ import annotations

import shutil
from importlib.resources import files
from pathlib import Path

from omegaconf import DictConfig, OmegaConf
from platformdirs import user_config_dir

_USER_CONFIG_DIR = Path(user_config_dir("pdbe_sifts"))
_USER_CONFIG_FILE = _USER_CONFIG_DIR / "config.yaml"


def get_template_path() -> Path:
    """Return the path to the built-in default config template."""
    return Path(
        files("pdbe_sifts.data").joinpath("default_config.yaml")
        # ._path
    )


def init_config(dest: Path | None = None) -> Path:
    """Copy the built-in config template to the user config directory.

    If dest is None, copies to ~/.config/pdbe_sifts/config.yaml
    """
    target = dest or _USER_CONFIG_FILE
    target.parent.mkdir(parents=True, exist_ok=True)

    if target.exists():
        print(f"Config already exists at: {target}")
        print("Use --force to overwrite.")
        return target

    shutil.copy(get_template_path(), target)
    print(f"Config template copied to: {target}")
    print("Please edit it before running the pipeline.")
    return target


def load_config(user_config: str | Path | None = None) -> DictConfig:
    """Load config with the following priority:
    1. Explicit --config path
    2. ~/.config/pdbe_sifts/config.yaml (set via pdbe_sifts init)
    3. Built-in defaults
    """
    # Built-in defaults
    cfg = OmegaConf.load(get_template_path())

    # Auto-detect user config if no explicit path given
    config_path = (
        Path(user_config)
        if user_config
        else (_USER_CONFIG_FILE if _USER_CONFIG_FILE.exists() else None)
    )

    if config_path is not None:
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        user_cfg = OmegaConf.load(config_path)
        cfg = OmegaConf.merge(cfg, user_cfg)

    return cfg
