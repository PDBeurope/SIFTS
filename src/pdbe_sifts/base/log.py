"""Logging for the whole package.

This module is used by all modules in the project for logging. It is
configured to output to the console. The log level is set to INFO by default,
but can be changed by setting the environment variable SIFTS_LOG_LEVEL to a
valid log level. The log level can be one of CRITICAL, ERROR, WARNING, INFO,
DEBUG, or NOTSET. NOTSET will disable logging.

Optionally uses coloredlogs package to color the output.

Example:
    >>> from pdbe_sifts.base.log import logger
    >>> logger.info("Hello world")
    >>> logger.error("Something went wrong", exc_info=True)
"""

import getpass
import logging
import os
import platform
import re

LOG_FORMAT = (
    "{hostname}: {username}: {asctime}: {module}: {lineno} {levelname}: {message}"
)

# Capture the original factory before overriding it.
old_factory = logging.getLogRecordFactory()


class SensitiveFormatter(logging.Formatter):
    """Formatter that removes sensitive information in URLs."""

    @staticmethod
    def _filter(s: str) -> str:
        return re.sub(r":\/\/(.*?)\@", r"://***:***/", s)

    def format(self, record: logging.LogRecord) -> str:  # noqa: A003
        original = logging.Formatter.format(self, record)
        return self._filter(original)


def record_factory(*args, **kwargs) -> logging.LogRecord:
    """Extend log records with hostname and username fields."""
    record = old_factory(*args, **kwargs)
    record.hostname = platform.node()
    record.username = getpass.getuser()
    return record


def _get_handler() -> logging.StreamHandler:
    """Build and return a console handler with the sensitive-data formatter.

    Returns:
        logging.StreamHandler: The handler to use for logging.
    """
    handler = logging.StreamHandler()
    handler.setFormatter(SensitiveFormatter(LOG_FORMAT, style="{"))
    return handler


logging.setLogRecordFactory(record_factory)

# Use a package-level logger so that third-party library messages are not
# captured, and propagation to the root logger is disabled to avoid duplicate
# output when coloredlogs is not installed.
logger = logging.getLogger("pdbe_sifts")
logger.propagate = False

log_level = os.environ.get("SIFTS_LOG_LEVEL", "INFO").upper()
try:
    logger.setLevel(log_level)
except ValueError:
    logger.error(f"Invalid log level '{log_level}' in SIFTS_LOG_LEVEL. Defaulting to INFO.")
    log_level = "INFO"
    logger.setLevel(log_level)

try:
    import coloredlogs

    # coloredlogs installs its own handler; we must NOT add ours as well or
    # every message would be printed twice.
    coloredlogs.install(level=log_level, logger=logger)
except ImportError:
    logger.addHandler(_get_handler())
