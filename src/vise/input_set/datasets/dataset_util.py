# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Dataset utilities for VASP input generation.

This module provides utilities for accessing VASP-related datasets and
performing calculations needed for input file generation.

Key Components:
    - LDAU: DFT+U parameter management
    - num_bands: Calculate number of bands for spectrum calculations
    - calc_kpar: Calculate optimal KPAR for k-point parallelization
    - incar_categories: Categorized INCAR flags from VASP
"""
from __future__ import annotations

import json
from itertools import chain
from math import ceil
from pathlib import Path
from typing import Any, Dict, List

from monty.serialization import loadfn
from pymatgen.core import Composition, Element
from pymatgen.io.vasp import Potcar, inputs

from vise.defaults import defaults

# Load INCAR parameters from pymatgen
_pymatgen_path = Path(inputs.__file__).parent
with open(_pymatgen_path / "incar_parameters.json", encoding="utf-8") as json_file:
    incar_params: Dict[str, Any] = json.loads(json_file.read())

#: Unoccupied band counts per element for NBANDS calculation
unoccupied_bands: Dict[str, int] = loadfn(Path(__file__).parent / "unoccupied_bands.yaml")

# Load INCAR flag categories
# Note: From Python 3.6+, dict maintains insertion order, so OrderedDict is not needed
incar_categories: Dict[str, List[str]] = dict(
    loadfn(Path(__file__).parent / "incar_flags.yaml")
)

# Add uncategorized tags to "others" category
_tag_set = set(chain.from_iterable(incar_categories.values()))
incar_categories["others"] = list(set(incar_params.keys()) - _tag_set)

# Add user-defined INCAR tag categories
for key, value in defaults.user_incar_tags.items():
    incar_categories[key] = value

#: Complete list of all INCAR flags across all categories
all_incar_flags: List[str] = sum(incar_categories.values(), [])


def has_f_elements(symbol_list: List[str]) -> bool:
    """Check if any element in the list has f-electrons (Z > 56).

    Args:
        symbol_list: List of element symbols.

    Returns:
        True if any element is a lanthanide or actinide (Z > 56).
    """
    return any([Element(el).Z > 56 for el in symbol_list])


class LDAU:
    """DFT+U parameter manager for Hubbard U corrections.

    This class loads U parameters from configuration files and provides
    the LDAUU, LDAUL, and LMAXMIX values needed for VASP DFT+U calculations.

    Attributes:
        symbol_list: List of element symbols in the structure.
        ldauu: List of U values for each element.
        ldaul: List of L values (orbital angular momentum) for each element.
    """

    def __init__(self, symbol_list: List[str]) -> None:
        """Initialize LDAU parameters from configuration.

        Args:
            symbol_list: List of element symbols in order of POSCAR.
        """
        self.symbol_list = symbol_list
        ldau_set = loadfn(defaults.u_parameter_set_yaml_file)
        self.ldauu: List[float] = [ldau_set["LDAUU"].get(el, 0) for el in symbol_list]
        self.ldaul: List[int] = [ldau_set["LDAUL"].get(el, -1) for el in symbol_list]

    @property
    def lmaxmix(self) -> int:
        """Determine LMAXMIX based on presence of f-elements.

        Returns:
            6 for f-element systems, 4 otherwise.
        """
        return 6 if has_f_elements(self.symbol_list) else 4

    @property
    def is_ldau_needed(self) -> bool:
        """Check if DFT+U corrections are needed.

        Returns:
            True if any element has a non-zero U value.
        """
        return set(self.ldauu) != {0}


def num_bands(
    composition: Composition,
    potcar: Potcar,
    spin_orbit: bool = False,
) -> int:
    """Calculate the number of bands for spectrum calculations.

    This function calculates NBANDS including unoccupied bands needed for
    optical absorption, band structure, and DOS calculations.

    Args:
        composition: pymatgen Composition object.
        potcar: pymatgen Potcar object.
        spin_orbit: Whether spin-orbit coupling is enabled (doubles bands).

    Returns:
        Total number of bands including unoccupied states.
    """
    total_bands = 0.0

    for element, potcar_single in zip(composition, potcar):
        num_atoms = composition[element]
        occupied_bands = potcar_single.nelectrons / 2
        bands_per_atom = occupied_bands + unoccupied_bands[str(element)]
        total_bands += num_atoms * bands_per_atom

    if spin_orbit:
        total_bands *= 2

    return ceil(total_bands)


def calc_kpar(
    num_kpoints: int,
    num_cores: int,
    unused_core_ratio_threshold: float,
) -> int:
    """Calculate optimal KPAR for k-point parallelization.

    KPAR determines how k-points are distributed across MPI processes.
    This function finds the largest valid KPAR that keeps the unused
    core ratio below the specified threshold.

    Args:
        num_kpoints: Total number of k-points in the calculation.
        num_cores: Number of available CPU cores.
        unused_core_ratio_threshold: Maximum acceptable ratio of unused cores.

    Returns:
        Optimal KPAR value (divisor of num_cores).

    Raises:
        ValueError: If no valid KPAR can satisfy the threshold constraint.
    """
    # Try divisors from largest to smallest
    divisors = [d for d in range(num_cores, 0, -1) if num_cores % d == 0]

    for divisor in divisors:
        kpt_per_core = round(num_kpoints / divisor, 5)
        fractional_part = kpt_per_core % 1

        # Calculate unused core ratio
        if fractional_part:
            fractional_part = 1 - fractional_part
        else:
            fractional_part = 0.0

        unused_core_ratio = fractional_part / ceil(kpt_per_core)

        if unused_core_ratio < unused_core_ratio_threshold:
            return divisor

    raise ValueError(
        f"The threshold for unused core ratio {unused_core_ratio_threshold} "
        f"is not adequate."
    )
