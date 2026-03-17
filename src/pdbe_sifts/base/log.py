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

LOG_FORMAT = "{hostname}: {username}: {asctime}: {module}: {lineno} {levelname}: {message}"


class SensitiveFormatter(logging.Formatter):
    """Formatter that removes sensitive information in urls."""

    @staticmethod
    def _filter(s):
        return re.sub(r":\/\/(.*?)\@", r"://***:***/", s)

    def format(self, record):  # noqa: A003
        original = logging.Formatter.format(self, record)
        return self._filter(original)


def record_factory(*args, **kwargs):
    record = old_factory(*args, **kwargs)
    record.hostname = platform.node()
    record.username = getpass.getuser()
    return record


def _get_handler() -> logging.StreamHandler:
    """Return a console (StreamHandler) for the package logger."""
    handler = logging.StreamHandler()
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

    coloredlogs.install(level=logLevel, logger=logger)
except ImportError:
    pass
