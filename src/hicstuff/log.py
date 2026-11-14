#!/usr/bin/env python3

"""Basic logging

Basic logging setup with three main handlers (stdout, log file and optionally
texting with the right API). By default the log file is disabled, but can be
enabled or changed using set_file_handler.
"""

import json
import logging
import logging.handlers

import requests

TEXT_CREDENTIALS_DEFAULT_PATH = "text_credentials.json"

CURRENT_LOG_LEVEL = logging.INFO


logging.captureWarnings(True)

logger = logging.getLogger()
logger.setLevel(CURRENT_LOG_LEVEL)

logfile_formatter = logging.Formatter(
    "%(asctime)s :: %(levelname)s :: %(message)s", datefmt="%Y-%m-%d,%H:%M:%S"
)
stdout_formatter = logging.Formatter("%(levelname)s :: %(message)s")
text_formatter = logging.Formatter("%(message)s")

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


def setup_text_logging(credentials=TEXT_CREDENTIALS_DEFAULT_PATH):
    """Setup text logging

    Setup text notifications on errors with a basic API. In order for it to
    work, a file named 'text_credentials.json' must be in the directory where
    a command is run. The logger performs a GET request to a text notification
    API.

    Parameters
    ----------
    credentials : str or pathlib.Path
        A JSON file with necessary credentials for the GET request
    """

    with open(credentials) as cred_handle:
        cred_dict = json.load(cred_handle)

    api_service = cred_dict.pop("api_service")

    class TextHandler(logging.Handler):
        def emit(self, record):
            log_entry = self.format(record)
            cred_dict["msg"] = log_entry
            return requests.get(api_service, cred_dict)

    sms_handler = TextHandler()
    sms_handler.setLevel(logging.ERROR)
    sms_handler.setFormatter(text_formatter)

    logger.addHandler(sms_handler)


try:
    setup_text_logging(TEXT_CREDENTIALS_DEFAULT_PATH)
except FileNotFoundError:
    pass
