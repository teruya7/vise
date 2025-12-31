# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Logging utilities for vise.

This module provides a configurable logger factory function that creates
loggers with consistent formatting across the vise package.
"""

import logging
import sys
from typing import IO, Literal, Optional

LoggerType = Literal["simple", "abundant"]


def get_logger(
    name: str,
    level: int = logging.DEBUG,
    stream: IO = sys.stdout,
    datetime_format: str = "%Y/%m/%d %H:%M:%S",
    logger_type: LoggerType = "simple",
    log_filename: Optional[str] = None,
    file_handler_level: int = logging.WARNING,
) -> logging.Logger:
    """Create and configure a logger with consistent formatting.

    Args:
        name: Logger name, typically __name__ of the calling module.
        level: Logging level for the logger (default: DEBUG).
        stream: Output stream for the stream handler (default: stdout).
        datetime_format: Format string for timestamps in log messages.
        logger_type: Format style - "simple" for minimal output,
                    "abundant" for detailed output with timestamps.
        log_filename: Optional path to a log file. If provided, a file
                     handler will be added.
        file_handler_level: Logging level for the file handler.

    Returns:
        Configured Logger instance.

    Raises:
        ValueError: If logger_type is not "simple" or "abundant".

    Examples:
        >>> logger = get_logger(__name__)
        >>> logger.info("Processing started")
    """
    # Define log formats
    if logger_type == "simple":
        log_format = "%(levelname)7s: %(message)s"
    elif logger_type == "abundant":
        log_format = "%(asctime)18s %(levelname)7s %(name)25s\n --> %(message)s"
    else:
        raise ValueError(
            f"Invalid logger_type: {logger_type!r}. "
            "Must be 'simple' or 'abundant'."
        )

    formatter = logging.Formatter(log_format, datefmt=datetime_format)
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Add file handler if filename is provided
    if log_filename:
        file_handler = logging.FileHandler(log_filename)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(file_handler_level)
        logger.addHandler(file_handler)

    # Add stream handler
    stream_handler = logging.StreamHandler(stream=stream)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger
