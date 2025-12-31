# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
"""POTCAR valence orbital utilities.

This module provides functions for extracting valence orbital information
from VASP POTCAR files.
"""

from typing import List

from pymatgen.io.vasp import PotcarSingle


def valence_orbitals_from_potcar(potcar_single: PotcarSingle) -> str:
    """Extract valence electron configuration from a POTCAR.

    Formats the electron configuration in a readable string format
    showing principal quantum number, orbital type, and electron count.

    Args:
        potcar_single: A pymatgen PotcarSingle object containing
                      the POTCAR data for a single element.

    Returns:
        Space-separated string of orbital configurations.
        Format: "(nl)m" where n=principal, l=orbital type, m=electrons.

    Examples:
        >>> # For an element with 3s2 3p6 configuration
        >>> valence_orbitals_from_potcar(potcar)
        '(3s)2 (3p)6'
    """
    orbital_strings: List[str] = []

    for config in potcar_single.electron_configuration:
        principal_num, orb_type, num_electrons = config
        orbital_strings.append(f"({principal_num}{orb_type}){num_electrons}")

    return " ".join(orbital_strings)
