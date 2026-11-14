#!/usr/bin/env python3

"""Basic logging

Basic logging setup with two handlers (stdout and log file). By default the
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


def set_file_handler(log_path, formatter=logfile_formatter):
    """Change the file handler for custom log file location"""

    filehandler = logging.FileHandler(log_path, "a")
    filehandler.setLevel(logging.INFO)
    filehandler.setFormatter(formatter)
    for hdlr in logger.handlers[:]:  # remove the existing file handlers
        if isinstance(hdlr, logging.FileHandler):
            logger.removeHandler(hdlr)
    logger.addHandler(filehandler)  # set the new handler
