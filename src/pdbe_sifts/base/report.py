"""Run reporting utilities (simplified logger-only version).

This is a standalone replacement for the ORC ``orc.base.report`` module.
All entry-level and summary-level events are written to the standard
:mod:`pdbe_sifts.base.log` logger; no database or SQLAlchemy dependency is
required.

Optionally, if the environment variable ``SIFTS_LOG_DB`` is set to a
SQLAlchemy-compatible connection string, structured logging to a database
can be enabled in a future extension of this module.
"""

from enum import Enum

from pdbe_sifts.base.log import logger


class EntryStatusType(Enum):
    """Status codes reported per processed entry or per batch run."""

    SUCCESS = "success"
    FAILED = "failed"
    SKIPPED = "skipped"
    KNOWN_EXCEPTION = "known_exception"

    def __repr__(self) -> str:
        return self.value


def log_entry(
    entry: str,
    module: str,
    status: EntryStatusType,
    exception_type=None,
    exception_message=None,
) -> None:
    """Log the processing result of a single entry.

    Args:
        entry: Entry identifier (e.g. PDB ID).
        module: Task class name.
        status: One of the :class:`EntryStatusType` values.
        exception_type: Exception class (``None`` on success).
        exception_message: Exception message string (``None`` on success).
    """
    if status == EntryStatusType.SUCCESS:
        logger.debug(f"[{module}] {entry} — {status}")
    else:
        logger.info(
            f"[{module}] {entry} — {status}"
            + (f" | {exception_type}: {exception_message}" if exception_type else "")
        )


def send_summary_message(
    module: str,
    status: EntryStatusType,
    processed: int | None = None,
    failed: int | None = None,
    known_exceptions: int | None = None,
    log_dir: str | None = None,
    message: str | None = None,
) -> None:
    """Log a summary for the completed batch run.

    Args:
        module: Task class name.
        status: Overall run status.
        processed: Total number of entries attempted.
        failed: Number of entries that failed (excluding known exceptions).
        known_exceptions: Number of entries skipped as known exceptions.
        log_dir: Directory where per-job logs were written.
        message: Optional free-form message.
    """
    parts = [
        f"[{module}] Batch complete — status={status}",
        f"processed={processed}",
        f"failed={failed}",
        f"known_exceptions={known_exceptions}",
    ]
    if log_dir:
        parts.append(f"log_dir={log_dir}")
    if message:
        parts.append(f"message={message}")

    summary = " | ".join(parts)
    if status == EntryStatusType.SUCCESS:
        logger.info(summary)
    else:
        logger.warning(summary)
