#!/usr/bin/env python3
from pathlib import Path
import os
import tomllib  # Python 3.11+ (use tomli for earlier versions)


class Config:
    """Global configuration for the PDBe-SIFTS pipeline."""

    def __init__(self, config_file: str | None = None):
        # Project root directory
        self.root = Path(__file__).resolve().parents[1]

        # Optional TOML file to override defaults
        self._config_data = {}
        config_path = Path(config_file) if config_file else self.root / "config.toml"
        if config_path.exists():
            with open(config_path, "rb") as f:
                self._config_data = tomllib.load(f)

        # Load subsections
        self.sifts = self._load_sifts_config()

    def _load_sifts_config(self):
        """Load paths and parameters related to SIFTS."""
        data = self._config_data.get("sifts", {})

        def resolve_path(value: str, default: str) -> Path:
            """Resolve a path: absolute paths are kept, relative ones are made relative to self.root."""
            p = Path(os.getenv(value, data.get(value.lower(), default)))
            return p if p.is_absolute() else (self.root / p)

        return type("SiftsConfig", (), {
            "tax_fix_pkl": resolve_path("SIFTS_TAX_FIX", "input_files/nf90_taxid.pkl"),
            "three_to_one": resolve_path(
                "SIFTS_THREE_TO_ONE",
                "input_files/three_to_one_letter_mapping.csv"
            ),
        })()


    def __repr__(self):
        return f"<Config root={self.root}>"
