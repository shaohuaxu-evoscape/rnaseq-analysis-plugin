"""Pipeline logging utilities."""

import logging
import sys


def setup_logger(name="ds002_pipeline", level=logging.INFO):
    """Create a formatted logger.

    Parameters
    ----------
    name : str
        Logger name.
    level : int
        Logging level.

    Returns
    -------
    logging.Logger
    """
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger

    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)-8s %(name)s | %(message)s",
        datefmt="%H:%M:%S",
    )
    handler.setFormatter(fmt)
    logger.addHandler(handler)
    logger.propagate = False
    return logger
