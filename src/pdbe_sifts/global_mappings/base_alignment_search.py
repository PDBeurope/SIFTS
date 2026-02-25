"""
Abstract base class for sequence alignment search tools.

This module defines the `AlignmentSearch` class, which provides a
standardized interface for running alignment searches between a query
and a target database. Subclasses must implement the `_process` method
to define the alignment logic using a specific tool (e.g., BLAST, MMseqs2, etc.).
"""

from abc import ABC, abstractmethod
from timeit import default_timer as timer
from pathlib import Path
from pdbe_sifts.base.log import logger

class AlignmentSearch(ABC):
    """
    Abstract base class for alignment search tools.

    This class defines a framework for running an alignment of a query
    against a target database and measuring its execution time. In subclasses,
    the _process method need to be implemented.

    Attributes:
        query_path (str): Path to the query file.
        target_path (str): Path to the target database or sequence file.
        output_path (str): Path where alignment results will be saved.
    """

    def __init__(self, query_path: str, target_path: str, output_path: str):
        """Initializes the AlignmentSearch with query, target, and output file paths."""
        self.query_path = query_path
        self.target_path = target_path
        self.output_path = output_path

        if not Path(self.query_path).exists():
            raise FileNotFoundError(f"Query file not found : {self.query_path}")

    def run(self):
        """Runs the alignment process and logs timing information."""
        logger.info(f"Processing the search of the query {self.query_path} against {self.target_path}.")
        start = timer()
        self._process()
        end = timer()
        logger.info("Search finished.")
        logger.info(f'Search process took: {end - start} seconds.')
        logger.info(f"Results saved: {self.output_path}")

    @abstractmethod
    def _process(self):
        """Defines how the alignment should be performed. Must be implemented by subclasses."""
        pass
