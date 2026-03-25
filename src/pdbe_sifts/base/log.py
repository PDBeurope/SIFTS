"""Logging for the whole package.

This module is used by all modules in the project for logging. It is
configured to output to the console. The log level is set to INFO by default,
but can be changed by setting the environment variable SIFTS_LOG_LEVEL to a
valid log level (CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET).

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
import sys

LOG_FORMAT = "{hostname}: {username}: {asctime}: {module}: {lineno} {levelname}: {message}"


class SensitiveFormatter(logging.Formatter):
    """Formatter that removes sensitive information in urls."""

    @staticmethod
    def _filter(s):
        """Remove credentials embedded in URLs of the form ``scheme://user:pass@host``.

        Args:
            s: Raw log string that may contain URLs with embedded credentials.

        Returns:
            The input string with any ``user:pass@`` segments replaced by
            ``***:***@``.
        """
        return re.sub(r":\/\/(.*?)\@", r"://***:***/", s)

    def format(self, record):  # noqa: A003
        """Format a log record and strip any embedded URL credentials.

        Args:
            record: The :class:`logging.LogRecord` to format.

        Returns:
            The fully formatted log string with credentials redacted.
        """
        original = logging.Formatter.format(self, record)
        return self._filter(original)


def record_factory(*args, **kwargs):
    """Extend the default :class:`logging.LogRecord` with hostname and username.

    Installed via :func:`logging.setLogRecordFactory` so that every log record
    automatically carries ``record.hostname`` and ``record.username`` attributes,
    making them available in format strings.

    Args:
        *args: Positional arguments forwarded to the original record factory.
        **kwargs: Keyword arguments forwarded to the original record factory.

    Returns:
        A :class:`logging.LogRecord` augmented with ``hostname`` and
        ``username`` fields.
    """
    record = old_factory(*args, **kwargs)
    record.hostname = platform.node()
    record.username = getpass.getuser()
    return record


def _get_handler() -> logging.StreamHandler:
    """Return a console (StreamHandler) for the package logger."""
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(SensitiveFormatter(LOG_FORMAT, style="{"))
    return handler


old_factory = logging.getLogRecordFactory()
logging.setLogRecordFactory(record_factory)

logger = logging.getLogger()
logger.addHandler(_get_handler())

envvar = "SIFTS_LOG_LEVEL"
logLevel = os.getenv(envvar, "INFO")

try:
    logger.setLevel(logLevel)
except ValueError:
    logger.error(f"Invalid log level in env var {envvar}. Defaulting to INFO")
    logLevel = "INFO"
    logger.setLevel(logLevel)

try:
    import coloredlogs

    coloredlogs.install(level=logLevel, logger=logger, stream=sys.stdout)
except ImportError:
    pass
