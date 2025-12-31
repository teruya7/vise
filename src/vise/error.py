# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Base exception classes for VISE.

This module defines the base exception hierarchy for the VISE package.
All VISE-specific exceptions should inherit from ViseError.
"""


class ViseError(Exception):
    """Base exception class for all VISE-related errors.

    All custom exceptions in VISE should inherit from this class to allow
    for unified exception handling across the package.

    Example:
        >>> raise ViseError("Something went wrong in VASP input generation")
    """

    pass
