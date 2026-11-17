#!/usr/bin/env python3

"""Basic logging

Basic logging setup with two handlers (standard error and log file). By default the
log file is disabled, but can be enabled or changed using set_file_handler.
"""

import logging
import logging.handlers

CURRENT_LOG_LEVEL = logging.INFO


logging.captureWarnings(True)

logger = logging.getLogger()
logger.setLevel(CURRENT_LOG_LEVEL)

logfile_formatter = logging.Formatter(
    "%(asctime)s :: %(levelname)s :: %(message)s", datefmt="%Y-%m-%d,%H:%M:%S"
)
stdout_formatter = logging.Formatter("%(levelname)s :: %(message)s")

# file_handler = logging.FileHandler("hicstuff.log", "a")

# file_handler.setLevel(logging.INFO)
# file_handler.setFormatter(logfile_formatter)
# logger.addHandler(file_handler)

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(stdout_formatter)
logger.addHandler(stream_handler)

# The log file handler set by set_file_handler, if any
_file_handler = None


def remove_file_handler():
    """Detach and close the log file handler set by set_file_handler, if
    any. Other handlers of the logger are left alone."""
    global _file_handler
    if _file_handler is not None:
        logger.removeHandler(_file_handler)
        _file_handler.close()
        _file_handler = None


def set_file_handler(log_path, formatter=logfile_formatter):
    """Change the file handler for custom log file location

    The handler set by a previous call is replaced. File handlers attached
    to the logger by other code (an application, pytest) are left alone.
    """
    global _file_handler
    filehandler = logging.FileHandler(log_path, "a")
    filehandler.setLevel(logging.INFO)
    filehandler.setFormatter(formatter)
    remove_file_handler()
    _file_handler = filehandler
    logger.addHandler(filehandler)  # set the new handler
