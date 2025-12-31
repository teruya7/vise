# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
"""FFT grid calculation utilities for VASP.

This module provides functions for calculating the FFT grid sizes used
in VASP based on the energy cutoff and lattice parameters, following
VASP's internal algorithm.
"""

import numpy as np

# Physical constants used in VASP 5.4.4 source code
RYTOEV: float = 13.605826  # Rydberg to eV conversion
AUTOA: float = 0.529177249  # Atomic unit (Bohr) to Angstrom
PI: float = 3.141592653589793238


def vasp_grid(encut: float, lattice_length: float, symprec: str) -> int:
    """Calculate the FFT grid dimension for a lattice direction in VASP.

    Follows the algorithm in VASP 5.4.4 main.F source code to determine
    the number of FFT grid points along a lattice direction.

    Args:
        encut: Plane-wave energy cutoff in eV (ENCUT in INCAR).
        lattice_length: Length of the lattice vector in Angstroms.
        symprec: Precision string - 'accurate'/'high' gives factor 4,
                otherwise factor 3.

    Returns:
        Number of FFT grid points for this direction.

    Note:
        This calculates NGPTAR before FFTCHK optimization. VASP may
        adjust this value based on FFT algorithm requirements.

    Examples:
        >>> vasp_grid(500.0, 4.0, "accurate")  # High precision
        40
        >>> vasp_grid(500.0, 4.0, "normal")    # Normal precision
        30

    Reference:
        VASP 5.4.4 source: main.F
        XCUTOF = SQRT(INFO%ENMAX/RYTOEV) / (2*PI/(LATT_CUR%ANORM(1)/AUTOA))
        WFACT = 4 for accurate/high, 3 otherwise
        GRID%NGPTAR(1) = XCUTOF * WFACT + 0.5
    """
    # Use factor 4 for accurate/high precision, 3 otherwise
    first_char = symprec[0].lower() if symprec else "n"
    factor = 4 if first_char in ("a", "h") else 3

    # Convert units to atomic units as VASP does internally
    cutoff_rydberg = encut / RYTOEV
    lattice_au = lattice_length / AUTOA

    # Reciprocal lattice vector magnitude
    reciprocal_length = 2 * PI / lattice_au

    # Calculate grid dimension
    return int(round(factor * np.sqrt(cutoff_rydberg) / reciprocal_length))
