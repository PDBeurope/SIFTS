from abc import ABC, abstractmethod
from datetime import datetime
from typing import Iterable
from enum import Enum

import gemmi

from pdbe_sifts.base.log import logger


class Modes(Enum):
    SINGLE = "SINGLE"
    BATCH = "BATCH"


class Batchable(ABC):
    """Abstract Base Class for tasks that can be parallelized."""

    used_cif_categories: Iterable[str] = set()

    @abstractmethod
    def process_entry(self, entry_id: str):
        """Abstract method to process single entry.

        Each subclass should implement this method to process a single entry.

        Args:
            entry_id: ID to process e.g. CCD ID, PDB ID, PRD ID, Uniprot Accession
        """
        pass


    def _is_future_date(self, date: str) -> bool:
        """Check if the date is in the future.

        Args:
            date (str): Date in the format "YYYY-MM-DD"

        Returns:
            bool: True if the date is in the future, False otherwise.
        """
        return datetime.strptime(date, "%Y-%m-%d") > datetime.now()

    def no_used_cif_category_modified(self, cif_file: str) -> bool:
        """Returns true if the CIF file has no modifications in the categories used by the task.

        Scenarios
        1. Missing _pdbx_audit_revision_history or _pdbx_audit_revision_category, return False
        2. `self.used_cif_categories` is empty, return False.
        3. Latest revision is not in the future, return False.
        4. Modified categories exist, but none of them are used by task, return True.

        Args:
            cif_file (str): Path to the CIF file

        Returns:
            bool: True if the CIF file has no modifications in the categories used by the task.
                  False otherwise including if the categories are not found in the CIF file.
        """
        if not self.used_cif_categories:
            logger.debug("No categories to check for modifications")
            return False

        # Remove leading _ in self.used_cif_categories if present
        # since the categories in pdbx_ausit_category are without leading _
        self.used_cif_categories = {cat.lstrip("_") for cat in self.used_cif_categories}

        block = gemmi.cif.read(str(cif_file)).sole_block()
        history = block.find(
            "_pdbx_audit_revision_history.", ["ordinal", "revision_date"]
        )
        ordinals = [
            row["ordinal"]
            for row in history
            if self._is_future_date(row["revision_date"])
        ]

        if not ordinals:
            logger.info("No future revisions found")
            return False

        categories = block.find(
            "_pdbx_audit_revision_category.", ["revision_ordinal", "category"]
        )
        if not categories:
            logger.info("No pdbx_audit_revision_category found")
            return False

        modified_categories = {
            row["category"] for row in categories if row["revision_ordinal"] in ordinals
        }
        modified_used_categories = modified_categories.intersection(
            self.used_cif_categories
        )
        if not modified_used_categories:
            logger.debug(f"Modified categories: {', '.join(modified_categories)}")
            logger.debug(f"Used categories: {', '.join(self.used_cif_categories)}")
            logger.info("None of the modified categories are used.")
            return True

        logger.info(f"Modified categories found: {', '.join(modified_used_categories)}")
        return False