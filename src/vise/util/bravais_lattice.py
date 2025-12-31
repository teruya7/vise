# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Bravais lattice types and utilities.

This module defines the BravaisLattice enum representing the 14 Bravais
lattice types and provides methods for determining lattice types from
space groups and k-point mesh configurations.
"""

from pathlib import Path
from typing import Dict, List, Tuple

from monty.serialization import loadfn

from vise.util.enum import ExtendedEnum

# Load Hermann-Mauguin symbols mapping: space group number -> symbol
_hermann_mauguin_yaml = Path(__file__).parent / "Hermann-Mauguin.yaml"
hermann_mauguin_list: Dict[int, str] = loadfn(_hermann_mauguin_yaml)

# Space group number ranges for each crystal system
CRYSTAL_SYSTEM_SG_RANGES: Dict[str, Tuple[int, int]] = {
    "a": (1, 2),      # Triclinic (anorthic)
    "m": (3, 15),     # Monoclinic
    "o": (16, 74),    # Orthorhombic
    "t": (75, 142),   # Tetragonal
    "h": (143, 194),  # Hexagonal/Trigonal
    "c": (195, 230),  # Cubic
}


class BravaisLattice(ExtendedEnum):
    """The 14 Bravais lattice types.

    Each lattice is named with a two-letter code:
    - First letter: crystal system (a, m, o, t, h, c)
    - Second letter: centering type (P, C, I, F, R, A)

    Attributes:
        aP: Triclinic primitive
        mP: Monoclinic primitive
        mC: Monoclinic C-centered
        oP: Orthorhombic primitive
        oF: Orthorhombic face-centered
        oI: Orthorhombic body-centered
        oA: Orthorhombic A-centered
        oC: Orthorhombic C-centered
        tP: Tetragonal primitive
        tI: Tetragonal body-centered
        hR: Hexagonal rhombohedral
        hP: Hexagonal primitive
        cP: Cubic primitive
        cF: Cubic face-centered
        cI: Cubic body-centered
    """

    aP = "aP"
    mP = "mP"
    mC = "mC"
    oP = "oP"
    oF = "oF"
    oI = "oI"
    oA = "oA"
    oC = "oC"
    tP = "tP"
    tI = "tI"
    hR = "hR"
    hP = "hP"
    cP = "cP"
    cF = "cF"
    cI = "cI"

    @classmethod
    def from_sg_num(cls, sg_num: int) -> "BravaisLattice":
        """Determine Bravais lattice type from space group number.

        Args:
            sg_num: Space group number (1-230).

        Returns:
            The corresponding BravaisLattice enum member.

        Raises:
            ValueError: If space group number is out of range.

        Examples:
            >>> BravaisLattice.from_sg_num(225)  # Fm-3m
            cF
        """
        # Get centering letter from Hermann-Mauguin symbol
        centering_letter = hermann_mauguin_list[sg_num][0]

        # Find crystal system from space group number
        for crystal_system, (sg_min, sg_max) in CRYSTAL_SYSTEM_SG_RANGES.items():
            if sg_min <= sg_num <= sg_max:
                bravais_code = crystal_system + centering_letter
                return cls.from_string(bravais_code)

        raise ValueError(f"Invalid space group number: {sg_num}")

    @property
    def kpt_centering(self) -> List[float]:
        """K-point mesh centering shift for this lattice type.

        Returns appropriate k-point shift for VASP calculations:
        - Gamma-centered [0, 0, 0] for oF, tI, cF
        - Half-shifted along c [0, 0, 0.5] for hP
        - Monkhorst-Pack [0.5, 0.5, 0.5] for others

        Returns:
            List of three floats representing the k-point shift.

        Note:
            cI does not require gamma centering despite being body-centered.
        """
        gamma_centered_lattices = [self.oF, self.tI, self.cF]
        if self in gamma_centered_lattices:
            return [0.0, 0.0, 0.0]
        elif self is self.hP:
            return [0.0, 0.0, 0.5]
        else:
            return [0.5, 0.5, 0.5]

    @property
    def need_same_num_kpt(self) -> bool:
        """Whether this lattice requires equal k-points in certain directions.

        Returns:
            True for oI and tI lattices, False otherwise.
        """
        return self in (self.oI, self.tI)
