# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
"""String formatting utilities for vise.

This module provides functions for formatting chemical formulas and
converting numbers to subscript characters.
"""

import re

# Mapping of digits to Unicode subscript characters
SUBSCRIPT_DIGITS: dict[int, str] = {
    0: "₀", 1: "₁", 2: "₂", 3: "₃", 4: "₄",
    5: "₅", 6: "₆", 7: "₇", 8: "₈", 9: "₉",
}


def latexify(formula: str) -> str:
    """Convert a chemical formula to LaTeX format with subscripts.

    Transforms element-number patterns (e.g., Fe2O3) to LaTeX subscript
    notation (e.g., Fe$_{2}$O$_{3}$).

    Note:
        This function was copied from pymatgen 2021.3.9 as it was
        scheduled for removal in 2022.

    Args:
        formula: Chemical formula string (e.g., "Fe2O3", "H2O").

    Returns:
        LaTeX-formatted formula with proper subscripts.

    Examples:
        >>> latexify("Fe2O3")
        'Fe$_{2}$O$_{3}$'
        >>> latexify("H2O")
        'H$_{2}$O'
    """
    return re.sub(r"([A-Za-z\(\)])(\d+\.?\d*)", r"\1$_{\2}$", formula)


def numbers_to_subscripts(text: str) -> str:
    """Convert all digits in a string to Unicode subscript characters.

    Args:
        text: Input string containing digits.

    Returns:
        String with all digits replaced by their subscript equivalents.

    Examples:
        >>> numbers_to_subscripts("H2O")
        'H₂O'
        >>> numbers_to_subscripts("Fe2O3")
        'Fe₂O₃'
    """
    result = text
    for digit, subscript in SUBSCRIPT_DIGITS.items():
        result = result.replace(str(digit), subscript)
    return result


# Backward compatibility alias
numbers_to_lowercases = numbers_to_subscripts
