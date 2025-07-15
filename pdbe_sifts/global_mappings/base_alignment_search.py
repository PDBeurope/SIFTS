from abc import ABC, abstractmethod
from timeit import default_timer as timer
from pathlib import Path
from pdbe_sifts.base.log import logger

class AlignmentSearch(ABC):
    def __init__(self, query_path: str, target_path: str, output_path: str):
        self.query_path = query_path
        self.target_path = target_path
        self.output_path = output_path

        if not Path(self.query_path).exists():
            raise FileNotFoundError(f"Query file not found : {self.query_path}")

    def run(self):
        """Run the alignment tool"""
        logger.info(f"Processing the search of the query {self.query_path} against {self.target_path}.")
        start = timer()
        self._process()
        end = timer()
        logger.info("Search finished.")
        logger.info(f'Search process took: {end - start} seconds.')
        logger.info(f"Results saved: {self.output_path}")

    @abstractmethod
    def _process(self):
        """Function to run the alignment tool. Implementation to be provided by subclasses"""
        pass
