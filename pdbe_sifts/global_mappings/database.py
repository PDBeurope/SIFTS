from abc import ABC, abstractmethod
from timeit import default_timer as timer
from pathlib import Path
from pdbe_sifts.base.log import logger

class ToolDatabase(ABC):
    def __init__(self, input_path: str, output_path: str):
        """"""
        self.input_path = Path(input_path)
        self.output_path = Path(output_path)

        if not self.input_path.exists():
            raise FileNotFoundError(f"Input file not found : {self.input_path}")

    def run(self):
        """Creates the database."""
        logger.info(f"Processing the creation of the database from {self.input_path}. This could take many hours/days.")
        start = timer()
        self._process()
        end = timer()
        logger.info("Creation finished.")
        logger.info(f'Creation process took: {end - start} seconds.')
        logger.info(f"Database saved: {self.output_path}")

    @abstractmethod
    def _process(self):
        """Function to create the database. Implementation to be provided by subclasses"""
        pass
