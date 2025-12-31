# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
"""Volumetric data handling for VASP.

This module provides utilities for processing volumetric data like
charge densities and wave functions from VASP calculations.
"""

from typing import List, Optional

import numpy as np
from numpy import prod
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Chgcar, Poscar, VolumetricData

from vise.util.logger import get_logger

logger = get_logger(__name__)

# Small value to avoid division by zero
_MINOR = 1e-3

# Default isosurface levels (fractions of maximum value)
DEFAULT_BORDER_FRACTIONS = [0.1, 0.5, 0.8]


def make_spin_charges(chgcar: Chgcar) -> List[Chgcar]:
    """Split spin-polarized CHGCAR into separate spin channels.

    Args:
        chgcar: CHGCAR with spin-polarized data.

    Returns:
        List of Chgcar objects, one per spin channel.
        For non-spin-polarized, returns single-element list.
    """
    result = [Chgcar(chgcar.structure, {"total": chgcar.spin_data[Spin.up]})]

    if "diff" in chgcar.data:
        result.append(
            Chgcar(chgcar.structure, {"total": chgcar.spin_data[Spin.down]})
        )

    return result


def light_weight_vol_text(
    volumetric_data: VolumetricData,
    border_fractions: Optional[List[float]] = None,
) -> str:
    """Create simplified volumetric data text representation.

    Converts continuous volumetric data to discrete levels based on
    threshold fractions, producing a more compact representation.

    Args:
        volumetric_data: VASP volumetric data (CHGCAR, LOCPOT, etc.).
        border_fractions: Threshold fractions of max value.
                         Points above each threshold increment by 1.

    Returns:
        Text representation suitable for lightweight visualization.
    """
    data = np.zeros(prod(volumetric_data.dim), dtype=int)

    # Normalize to [0, 1] with small offset to avoid zero
    normalized = (
        volumetric_data.data["total"] / np.max(volumetric_data.data["total"])
    ) + _MINOR

    border_fractions = border_fractions or DEFAULT_BORDER_FRACTIONS

    for border in border_fractions:
        # Transpose needed as VASP uses column-major (Fortran) order
        data += (normalized > border).T.flatten()

    lines = [
        Poscar(volumetric_data.structure).get_str(),
        " ".join(str(d) for d in volumetric_data.dim),
        " ".join(data.astype(str)),
    ]

    return "\n".join(lines)


# Backward compatibility aliases
default_border_fractions = DEFAULT_BORDER_FRACTIONS
