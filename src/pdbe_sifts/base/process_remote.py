#!/usr/bin/env python3
import importlib
import pickle  # nosec
import sys


class _PipelineUnpickler(pickle.Unpickler):
    """Resolves __main__.ClassName to the actual pipeline module.

    When a Batchable subclass is pickled from a script run as __main__,
    pickle stores the class module as '__main__'. This unpickler redirects
    those lookups to the real module at unpickle time (no circular import).
    """
    _SEARCH_MODULES = [
        'pdbe_sifts.sifts_segments_generation',
        # Add other pipeline modules here if needed
    ]

    def find_class(self, module, name):
        if module == '__main__':
            for mod_name in self._SEARCH_MODULES:
                try:
                    mod = importlib.import_module(mod_name)
                    if hasattr(mod, name):
                        return getattr(mod, name)
                except ImportError:
                    pass
        return super().find_class(module, name)


def run():
    """Executes the process_remote method of the pickled object passed in."""
    if len(sys.argv) != 2:
        raise Exception(f"Invalid usage: \n{sys.argv[0]} </path/to/pickle_file>")
    obj = load_pickle(sys.argv[1])
    obj.process_remote()


def load_pickle(file_name: str) -> object:
    """Loads a pickle file and returns the object.

    Args:
        file_name (str): Path to pickle file to load

    Returns:
        Pickled object
    """
    with open(file_name, "rb") as p:
        return _PipelineUnpickler(p).load()  # nosec


if __name__ == "__main__":
    run()
