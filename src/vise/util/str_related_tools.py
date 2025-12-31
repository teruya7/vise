# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""String-to-type conversion utilities.

This module provides helper functions for converting strings to
boolean values and checking if strings represent numeric types.
"""

from vise.util.logger import get_logger

logger = get_logger(__name__)


def str2bool(value: str) -> bool:
    """Convert a string representation to a boolean value.

    Args:
        value: String to convert. Accepts 't', 'true' (case-insensitive)
               for True, and 'false' for False.

    Returns:
        Boolean value corresponding to the input string.

    Raises:
        ValueError: If the string is not a recognized boolean representation.

    Examples:
        >>> str2bool("true")
        True
        >>> str2bool("FALSE")
        False
        >>> str2bool("T")
        True
    """
    normalized = value.lower()
    if normalized in ("t", "true"):
        return True
    elif normalized in ("false",):
        return False
    else:
        raise ValueError("invalid truth value %r" % (value,))



def is_str_digit(text: str) -> bool:
    """Check if a string represents a valid number (int or float).

    Args:
        text: String to check.

    Returns:
        True if the string can be parsed as a float, False otherwise.

    Examples:
        >>> is_str_digit("123")
        True
        >>> is_str_digit("3.14")
        True
        >>> is_str_digit("abc")
        False
    """
    try:
        float(text)
        return True
    except ValueError:
        return False


def is_str_int(text: str, rounding_error: float = 1e-7) -> bool:
    """Check if a string represents an integer value.

    Args:
        text: String to check.
        rounding_error: Maximum allowed difference from an integer value.

    Returns:
        True if the string represents an integer within tolerance,
        False otherwise.

    Examples:
        >>> is_str_int("42")
        True
        >>> is_str_int("3.14")
        False
    """
    try:
        if int(text) - float(text) < rounding_error:
            return True
        else:
            return False
    except ValueError:
        return False

