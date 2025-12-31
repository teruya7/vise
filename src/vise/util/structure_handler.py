# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Structure handling utilities for vise.

This module provides helper functions for manipulating crystal structures,
including symbol extraction and matrix sanitization.
"""

from itertools import groupby
from typing import List, Union

import numpy as np
from pymatgen.core import Structure

from vise.util.logger import get_logger

logger = get_logger(__name__)


def create_symbol_list(structure: Structure) -> List[str]:
    """Extract unique consecutive element symbols from structure.

    Removes consecutive duplicates while preserving order and handling
    non-consecutive repetitions.

    Args:
        structure: pymatgen Structure object.

    Returns:
        List of element symbols with consecutive duplicates removed.

    Examples:
        >>> # For structure with species ["H", "H", "O", "O", "H"]
        >>> create_symbol_list(structure)
        ['H', 'O', 'H']
    """
    species = [str(s) for s in structure.species]
    # Use itertools groupby recipe: unique_justseen
    # https://docs.python.org/3/library/itertools.html#itertools-recipes
    return [key for key, _ in groupby(species)]


def sanitize_matrix(
    matrix: Union[List[int], List[List[int]]]
) -> List[List[int]]:
    """Convert various matrix input formats to a 3x3 matrix.

    Supports several input formats:
    - 9 elements: Reshape to 3x3 matrix
    - 3 elements: Create diagonal matrix
    - 1 element: Create uniform diagonal matrix

    Args:
        matrix: Matrix elements in one of the supported formats.

    Returns:
        3x3 matrix as nested list.

    Raises:
        ValueError: If matrix has unsupported number of elements.

    Examples:
        >>> sanitize_matrix([2, 2, 2])
        [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        >>> sanitize_matrix([3])
        [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
    """
    length = len(matrix)

    if length == 9:
        # Reshape flat list to 3x3
        return [matrix[:3], matrix[3:6], matrix[6:]]
    elif length == 3:
        # Create diagonal matrix from 3 values
        result = np.eye(3, dtype=int)
        for i in range(3):
            result[i, i] = matrix[i]
        return result.tolist()
    elif length == 1:
        # Create uniform diagonal matrix
        result = np.eye(3, dtype=int) * matrix[0]
        return result.tolist()
    else:
        raise ValueError(
            f"Unsupported matrix length: {length}. "
            "Expected 1, 3, or 9 elements."
        )