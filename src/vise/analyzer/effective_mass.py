# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Effective mass calculation and analysis.

This module provides the EffectiveMass dataclass for storing and analyzing
carrier effective masses from transport calculations (typically BoltzTraP2).
"""

from dataclasses import dataclass
from typing import List, Tuple, Union

import numpy as np
from monty.json import MSONable
from numpy import linalg
from tabulate import tabulate

from vise.util.mix_in import ToJsonFileMixIn

# Type aliases
Matrix3x3 = List[List[float]]


@dataclass
class EffectiveMass(MSONable, ToJsonFileMixIn):
    """Carrier effective mass data from transport calculations.

    Stores effective mass tensors for holes (p) and electrons (n) at
    different carrier concentrations and a fixed temperature.

    Attributes:
        p: Hole effective mass tensors [concentration][3x3 matrix].
        n: Electron effective mass tensors [concentration][3x3 matrix].
        temperature: Temperature for the calculation (K).
        concentrations: List of carrier concentrations (cm^-3).
    """

    p: List[Matrix3x3]
    n: List[Matrix3x3]
    temperature: float
    concentrations: List[float]

    def average_mass(
        self, carrier_type: str, concentration: float
    ) -> float:
        """Calculate average effective mass (trace / 3).

        Args:
            carrier_type: 'p' for holes, 'n' for electrons.
            concentration: Carrier concentration to query.

        Returns:
            Average of diagonal elements of mass tensor.
        """
        mass_tensor = self.effective_mass(carrier_type, concentration)
        return sum(mass_tensor[i][i] for i in range(3)) / 3

    def minimum_mass(
        self, carrier_type: str, concentration: float
    ) -> float:
        """Calculate minimum effective mass (smallest eigenvalue).

        Args:
            carrier_type: 'p' for holes, 'n' for electrons.
            concentration: Carrier concentration to query.

        Returns:
            Minimum eigenvalue of the mass tensor.

        Raises:
            ValueError: If minimum mass has significant imaginary component.
        """
        mass_tensor = self.effective_mass(carrier_type, concentration)
        min_eigval, _ = lowest_eigval_and_vecs(np.array(mass_tensor))

        if isinstance(min_eigval, complex):
            if abs(min_eigval.imag) < 1e-3:
                return min_eigval.real
            else:
                raise ValueError(
                    f"Minimum effective mass has complex value: {min_eigval}"
                )
        return min_eigval

    def effective_mass(
        self, carrier_type: str, concentration: float
    ) -> Matrix3x3:
        """Get effective mass tensor for specified carrier and concentration.

        Args:
            carrier_type: 'p' for holes, 'n' for electrons.
            concentration: Carrier concentration to query.

        Returns:
            3x3 effective mass tensor in units of electron mass.
        """
        conc_index = self.concentrations.index(concentration)
        return getattr(self, carrier_type)[conc_index]

    def __str__(self) -> str:
        """Format effective masses for display."""
        lines = [f"temperature: {self.temperature}"]

        for conc, hole_mass, elec_mass in zip(
            self.concentrations, self.p, self.n
        ):
            lines.append("-" * 30)
            lines.append(f"concentration: {conc:g}")
            lines.append("p:")
            lines.append(tabulate(hole_mass))
            lines.append("n:")
            lines.append(tabulate(elec_mass))

        return "\n".join(lines)


def eigvals_and_vecs(
    matrix: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate eigenvalues and eigenvectors, sorted by eigenvalue.

    Args:
        matrix: 3x3 matrix to diagonalize.

    Returns:
        Tuple of (sorted eigenvalues, corresponding eigenvectors).
    """
    assert matrix.shape == (3, 3), "Matrix must be 3x3"

    eigvals, eigvecs = linalg.eig(matrix)
    sort_idx = eigvals.argsort()
    eigvals = eigvals[sort_idx]
    eigvecs = eigvecs[:, sort_idx].T  # Transpose for row-major format

    return eigvals, eigvecs


def lowest_eigval_and_vecs(
    matrix: np.ndarray,
) -> Tuple[Union[float, complex], np.ndarray]:
    """Find lowest eigenvalue and associated eigenvectors.

    Finds all eigenvectors with eigenvalues within 0.01% of the minimum.

    Args:
        matrix: 3x3 matrix to analyze.

    Returns:
        Tuple of (lowest eigenvalue, eigenvectors with that eigenvalue).
    """
    eigvals, eigvecs = eigvals_and_vecs(matrix)

    assert eigvals[0] > 0.0, "Lowest eigenvalue must be positive"

    min_eigval = eigvals[0]
    # Find all eigenvalues close to minimum (within 0.01%)
    near_min_mask = eigvals < min_eigval * 1.0001
    associated_eigvecs = np.round(eigvecs[near_min_mask], 2)

    return min_eigval, associated_eigvecs
