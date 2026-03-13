"""Logging for the whole package.

This module is used by all modules in the project for logging. It is
configured to output to the console. The log level is set to INFO by default, but can
be changed by setting the environment variable SIFTS_LOG_LEVEL to a valid
log level. The log level can be one of CRITICAL, ERROR, WARNING, INFO,
DEBUG, or NOTSET. NOTSET will disable logging.

Optionally uses coloredlogs package to color the output.

Example:
    >>> from pdbe_sifts.base.log import logger
    >>> logger.info("Hello world")
    >>> logger.error("Something went wrong", exc_info=True)
"""

import getpass
import logging
import logging.config
import os
import platform
import re

from pdbe_sifts.base.paths import (
    get_conf_logging_destination,
    get_conf_logging_elastic_file,
    get_conf_logging_file_path,
    get_conf_logging_max_bytes,
)

LOG_FORMAT = (
    "{hostname}: {username}: {asctime}: {module}: {lineno} {levelname}: {message}"
)


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


def _get_handler():
    """Get the handler used by the logger.

    The handler is determined by the environment variable SIFTS_LOG_DESTINATION.
    If SIFTS_LOG_DESTINATION is not set, the value in the conf.logging.destination
    is used.

    The possible values are:
    - console: Log to the console
    - file: Log to a file
    - syslog: Log to syslog
    - elastic: Log to elastic

    Returns:
        logging.Handler: The handler to use for logging.
    """

    log_destination = os.getenv("SIFTS_LOG_DESTINATION", get_conf_logging_destination())

    if log_destination == "file":
        handler = logging.RotatingFileHandler(
            get_conf_logging_file_path(), backupCount=10, maxBytes=get_conf_logging_max_bytes()
        )
        handler.setFormatter(SensitiveFormatter(LOG_FORMAT, style="{"))
    elif log_destination == "syslog":
        handler = logging.handlers.SysLogHandler()
        handler.setFormatter(SensitiveFormatter(LOG_FORMAT, style="{"))
    elif log_destination == "elastic":
        handler = logging.FileHandler(get_conf_logging_elastic_file())
        handler.setFormatter(SensitiveFormatter(LOG_FORMAT, style="{"))

    else:
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
