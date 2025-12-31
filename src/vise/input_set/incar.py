# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""INCAR file handling with pretty formatting.

This module provides ViseIncar, an enhanced version of pymatgen's Incar
class with improved parsing and formatted output grouped by category.
"""

import os
import re
from collections import OrderedDict
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List

from monty.io import zopen
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.vasp import Incar
from tabulate import tabulate

from vise.input_set.datasets.dataset_util import incar_categories
from vise.util.logger import get_logger

MODULE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))

logger = get_logger(__name__)


class ViseIncar(Incar):
    """Enhanced INCAR class with categorized formatting.

    Extends pymatgen's Incar to provide:
    - Pretty-printed output grouped by category
    - Support for ";" separator in INCAR lines
    - Improved parsing and serialization
    """

    @classmethod
    def from_file(cls, filename: str) -> "ViseIncar":
        """Read ViseIncar from file.

        Args:
            filename: Path to INCAR file.

        Returns:
            Parsed ViseIncar object.
        """
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read())

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ViseIncar":
        """Create ViseIncar from dictionary.

        Args:
            d: Dictionary of INCAR parameters.

        Returns:
            ViseIncar object.
        """
        kwargs = deepcopy(d)

        # Handle MAGMOM as list of Magmom objects
        if kwargs.get("MAGMOM") and isinstance(kwargs["MAGMOM"][0], dict):
            kwargs["MAGMOM"] = [Magmom.from_dict(m) for m in kwargs["MAGMOM"]]

        return cls({
            k: v for k, v in d.items()
            if k not in ("@module", "@class")
        })

    @classmethod
    def from_string(cls, string: str) -> "ViseIncar":
        """Parse INCAR from string content.

        Supports ";" separator for multiple parameters on one line,
        which is valid VASP syntax but not handled by base pymatgen.

        Args:
            string: INCAR file content.

        Returns:
            Parsed ViseIncar object.
        """
        params: Dict[str, Any] = {}

        for line in string.splitlines():
            # Support ";" semantic for splitting flags
            for sub_line in line.split(";"):
                match = re.match(r"\s*(\w+)\s*=\s*(.*)\s*", sub_line)
                if match:
                    key = match.group(1).strip()
                    val = match.group(2).strip()
                    params[key] = ViseIncar.proc_val(key, val)

        return cls(params)

    def __add__(self, other: Incar) -> "ViseIncar":
        """Merge two INCAR objects.

        Args:
            other: Another Incar to merge.

        Returns:
            New ViseIncar with combined settings.

        Raises:
            ValueError: If INCARs have conflicting values.
        """
        merged = dict(self.items())

        for key, value in other.items():
            if key in self and value != self[key]:
                raise ValueError(
                    f"INCARs have conflicting values for {key}: "
                    f"{self[key]} vs {value}"
                )
            merged[key] = value

        return ViseIncar(merged)

    def get_string(self, **kwargs) -> str:
        """Generate formatted INCAR string grouped by category.

        Returns:
            INCAR content with parameters organized by functional category.
        """
        sections: Dict[str, str] = OrderedDict()
        remaining_flags = list(self.keys())

        # Group parameters by category
        for category, category_flags in incar_categories.items():
            category_settings: List[List[str]] = []

            for flag in category_flags:
                if flag in remaining_flags:
                    category_settings.append([flag, self._setting_to_str(flag)])
                    remaining_flags.remove(flag)

            if category_settings:
                sections[f"# {category}"] = _tabulated_string(category_settings)

        # Warn about unrecognized flags
        if remaining_flags:
            logger.error(
                f"Unknown INCAR flags: {remaining_flags}. "
                f"These may be invalid; please add manually if needed."
            )

        return "\n\n".join(
            "\n".join([header, content])
            for header, content in sections.items()
        )

    def get_str(self, **kwargs) -> str:
        """Alias for get_string for compatibility."""
        return self.get_string(**kwargs)

    def _setting_to_str(self, flag: str) -> str:
        """Convert INCAR setting to string representation."""
        value = self[flag]
        if isinstance(value, list):
            return " ".join(str(item) for item in value)
        return str(value)

    # Keep old name for backward compatibility
    setting_to_str = _setting_to_str

    @property
    def is_ncl_calc(self) -> bool:
        """Check if calculation is non-collinear."""
        return any([
            self.get("LNONCOLLINEAR", False),
            self.get("LSORBIT", False),
        ])


def _tabulated_string(settings: List[List[str]]) -> str:
    """Format key-value pairs as aligned table with "=" separator.

    Args:
        settings: List of [key, value] pairs.

    Returns:
        Formatted string with aligned columns.

    Examples:
        >>> _tabulated_string([["NSW", "1"], ["ISMEAR", "-5"]])
        'NSW    = 1\\nISMEAR = -5'
    """
    # Add "=" separator column
    rows = [[key, "=", value] for key, value in settings]
    return str(tabulate(rows, tablefmt="plain", disable_numparse=True))


# Backward compatibility alias
tabulated_string = _tabulated_string
