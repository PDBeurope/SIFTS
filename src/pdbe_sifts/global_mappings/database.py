
"""
Base class for tool-specific database creation.

This module defines the abstract class `ToolDatabase`, which provides a
common interface and workflow for building structured databases from input data.
Subclasses must implement the `_process` method.
"""

from abc import ABC, abstractmethod
from timeit import default_timer as timer
from pathlib import Path
from pdbe_sifts.base.log import logger

class ToolDatabase(ABC):
    """
    Abstract base class for creating and managing tool-specific databases.

    Subclasses must implement the `_process` method, which defines how
    the database is built from the given input file. The class handles
    input validation, timing, and logging of the creation process.

    Attributes:
        input_path (Path): Path to the input file used to create the database.
        output_path (Path): Path where the resulting database will be saved.
    """

    def __init__(self, input_path: str, output_path: str):
        """Initializes the ToolDatabase with input and output file paths."""
        self.input_path = Path(input_path)
        self.output_path = Path(output_path)

        if not self.input_path.exists():
            raise FileNotFoundError(f"Input file not found : {self.input_path}")

    def run(self):
        """Creates the database and logs the processing time."""
        logger.info(f"Processing the creation of the database from {self.input_path}. This could take many hours/days.")
        start = timer()
        self._process()
        end = timer()
        logger.info("Creation finished.")
        logger.info(f'Creation process took: {end - start} seconds.')
        logger.info(f"Database saved: {self.output_path}")

    @abstractmethod
    def _process(self):
        """Defines how the database should be created. Must be implemented by subclasses."""
        pass
