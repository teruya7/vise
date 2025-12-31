# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
"""Irreducible representation extraction from VASP.

This module provides utilities for extracting irreducible representations
(irreps) from VASP band structure calculations using the irrep package.
"""

from typing import List, Tuple

from pymatgen.io.vasp import Kpoints

from vise.analyzer.plot_band import Irrep, Irreps
from vise.analyzer.vasp.plot_band import greek_to_unicode
from vise.error import ViseError
from vise.util.logger import get_logger

logger = get_logger(__name__)


def special_points_from_kpoints(
    kpoints_filename: str,
) -> Tuple[List[str], List[int]]:
    """Extract special point labels and indices from KPOINTS file.

    Args:
        kpoints_filename: Path to VASP KPOINTS file.

    Returns:
        Tuple of (special_point_labels, k-point_indices).
        Labels are converted from VASP convention (e.g., "GAMMA" -> "GM").
    """
    kpoints = Kpoints.from_file(kpoints_filename)
    special_points: List[str] = []
    kpt_indices: List[int] = []

    for idx, label in enumerate(kpoints.labels, 1):
        if label in [None, "None"] or label in special_points:
            continue
        special_points.append(label)
        kpt_indices.append(idx)

    # Convert GAMMA to GM for irrep package compatibility
    special_points = [x.replace("GAMMA", "GM") for x in special_points]

    return special_points, kpt_indices


# NOTE: make_irreps_from_wavecar is commented out as it requires
# the external 'irrep' package which may not be installed.
# See original commented code for the full implementation.


class ViseNoIrrepError(ViseError):
    """Error raised when irrep cannot be determined."""
    pass


def find_irrep(d: dict, threshold: float = 0.99) -> str:
    """Find the dominant irreducible representation.

    Args:
        d: Dictionary mapping irrep symbols to (weight, ...) tuples.
        threshold: Minimum weight to accept an irrep.

    Returns:
        Symbol of the dominant irrep.

    Raises:
        ViseNoIrrepError: If no irrep exceeds the threshold.
    """
    for symbol, values in d.items():
        if values[0] > threshold:
            return symbol

    logger.warning(
        f"No irrep found above threshold.\n"
        f"Threshold: {threshold}\n"
        f"Characters: {d}"
    )
    raise ViseNoIrrepError
