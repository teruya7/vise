# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""POTCAR set management for different calculation types.

This module provides the PotcarSet enum which defines different POTCAR
selection presets for various calculation types (normal, Materials Project,
GW calculations).
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Dict, Optional

from monty.serialization import loadfn

from vise.util.enum import ExtendedEnum


class PotcarSet(ExtendedEnum):
    """Enum defining POTCAR selection presets for different calculation types.

    Each preset maps element symbols to specific POTCAR variants based on
    the calculation requirements:
    - normal: Standard calculations with balanced accuracy/speed
    - mp_relax_set: Materials Project compatibility
    - gw: GW calculations requiring more accurate potentials

    Attributes:
        normal: Standard POTCAR selection.
        mp_relax_set: Materials Project compatible selection.
        gw: GW-appropriate POTCAR selection.
    """

    normal = "normal"
    mp_relax_set = "mp_relax_set"
    gw = "gw"

    def overridden_potcar_dict(
        self,
        override_potcar_set: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """Get POTCAR mapping with optional overrides.

        Args:
            override_potcar_set: Optional dictionary of element-to-POTCAR
                mappings to override the default selection.

        Returns:
            Dictionary mapping element symbols to POTCAR names.
            Example: {"Zr": "Zr_pv", "O": "O", ...}
        """
        result = deepcopy(self.potcar_dict())
        if override_potcar_set:
            result.update(override_potcar_set)
        return result

    def potcar_dict(self) -> Dict[str, str]:
        """Get the POTCAR mapping for this preset.

        Loads POTCAR selections from the potcar_set.yaml configuration file
        and returns the mapping for the current preset.

        Returns:
            Dictionary mapping element symbols to POTCAR names.
        """
        potcar_list = loadfn(Path(__file__).parent / "potcar_set.yaml")
        set_names = potcar_list.pop("set_names")

        # Initialize result dictionaries for each preset
        result: Dict[str, Dict[str, Optional[str]]] = {}
        for name in set_names:
            result[name] = {}

        # Parse POTCAR selections for each element
        for element, potcar_string in potcar_list.items():
            potcars = potcar_string.split()

            def sanitize(val: str) -> Optional[str]:
                """Convert special values to appropriate Python types.

                Args:
                    val: Raw value from YAML file.

                Returns:
                    Sanitized POTCAR name or None.
                """
                if val == "---":
                    # Use the first (default) POTCAR
                    return potcars[0]
                elif val == "None":
                    return None
                else:
                    return val

            for index, potcar_single in enumerate(potcars):
                set_name = set_names[index]
                result[set_name][element] = sanitize(potcar_single)

        return result[self.value]
