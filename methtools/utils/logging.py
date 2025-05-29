import logging
import sys


def setup_logging(level=logging.INFO):
    """
    Configures the logging to write messages to stderr with the specified level.

    Parameters:
        level (int): Logging level (e.g., logging.INFO, logging.DEBUG).
    """
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(level)  # Set the logging level dynamically

    # Clear existing handlers to avoid duplicate logs
    if logger.hasHandlers():
        logger.handlers.clear()

    # Define a logging format
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # Create a StreamHandler for stderr
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(level)  # Set handler level
    stderr_handler.setFormatter(formatter)  # Set formatter

    # Add the handler to the logger
    logger.addHandler(stderr_handler)
