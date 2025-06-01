import logging
import sys
from typing import Optional

def setup_logging(level: int = logging.INFO, log_file: Optional[str] = None) -> None:
    """
    Configures logging to write messages to stderr and optionally to a file.

    Parameters:
        level (int): The logging level to set for the root logger 
                     (e.g., logging.INFO, logging.DEBUG). Defaults to logging.INFO.
        log_file (Optional[str]): If provided, logs will also be written to this
                                  file path. Defaults to None.
    """
    # Get the root logger.
    # Subsequent calls to logging.getLogger('some_module') will inherit this config.
    logger = logging.getLogger()
    logger.setLevel(level)

    # Clear any existing handlers to avoid duplicate log messages in case
    # this function is called multiple times.
    if logger.hasHandlers():
        logger.handlers.clear()

    # Define a more informative logging format.
    # %(name)s will show the logger's name (e.g., 'methtools.mhb.caller')
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # --- Create Handler for Console (stderr) ---
    # This handler is always added.
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(level)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    # --- Create Handler for File (optional) ---
    # This handler is only added if a log_file path is provided.
    if log_file:
        try:
            # Use 'a' for append mode, so logs aren't erased on each run.
            file_handler = logging.FileHandler(log_file, mode='a')
            file_handler.setLevel(level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            # Log that we've successfully set up file logging.
            # This message will go to both console and the file.
            logging.info(f"Logging configured to also write to file: {log_file}")
        except (IOError, OSError) as e:
            # If the file can't be opened, log an error to the console and continue.
            logging.error(f"Failed to set up log file at '{log_file}': {e}", exc_info=True)


# --- Example Usage ---
# You would call this from your main cli.py

if __name__ == '__main__':
    # --- Backward-compatible call (works exactly as before) ---
    print("--- Running with default logging (console only) ---")
    setup_logging(level=logging.DEBUG)
    logging.info("This is an info message.")
    logging.debug("This is a debug message.")

    # --- New enhanced call (logs to console AND file) ---
    print("\n--- Running with file logging enabled ---")
    setup_logging(level=logging.INFO, log_file="my_app.log")
    logging.info("This message will go to the console and my_app.log.")
    logging.warning("This is a warning.")
    
    # Example of a logger in a submodule
    sub_logger = logging.getLogger("methtools.sub.module")
    sub_logger.info("This log message will show its specific origin.")

    print(f"\nCheck the contents of 'my_app.log' to see the output.")