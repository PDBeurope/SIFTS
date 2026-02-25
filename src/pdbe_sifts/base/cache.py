import os
from pathlib import Path
from platformdirs import user_cache_dir

_CONFIG_FILE = Path.home() / ".config" / "pdbe_sifts" / "unp_cache_dir.txt"


def set_unp_cache(path: Path) -> None:
    """Persist the UNP cache path to the user config file.
    Called once via: pdbe_sifts set-unp-cache /path/to/dir
    """
    _CONFIG_FILE.parent.mkdir(parents=True, exist_ok=True)
    _CONFIG_FILE.write_text(str(path.resolve()))
    print(f"UNP cache set to: {path.resolve()}")


def _get_unp_cache_path() -> Path:
    """Return the UNP cache directory path.

    Priority:
    1. Environment variable UNP_CACHE
    2. Persistent config file (~/.config/pdbe_sifts/unp_cache_dir.txt)
    3. Default platformdirs (~/.cache/pdbe_sifts/)
    """
    env_dir = os.environ.get("UNP_CACHE")
    if env_dir:
        unp_cache_dir = Path(env_dir)
    elif _CONFIG_FILE.exists():
        unp_cache_dir = Path(_CONFIG_FILE.read_text().strip())
    else:
        unp_cache_dir = Path(user_cache_dir("pdbe_sifts"))

    unp_cache_dir.mkdir(parents=True, exist_ok=True)
    return unp_cache_dir
