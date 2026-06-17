from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("pdbe_sifts")
except PackageNotFoundError:
    __version__ = "1.0"

__all__ = ["__version__"]
