#!/usr/bin/env python3


import pickle  # nosec
import sys


def run():
    """Executes the process_redis method of the pickled object passed in."""
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
        return pickle.load(p)  # nosec


if __name__ == "__main__":
    run()
