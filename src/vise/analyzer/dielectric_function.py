# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Dielectric function analysis and optical properties.

This module provides classes and functions for processing dielectric
function data and deriving optical properties like absorption coefficient,
refractive index, and reflectivity.
"""

from dataclasses import dataclass
from math import pi, sqrt
from typing import List, Optional

import numpy as np
import pandas as pd
from monty.json import MSONable
from scipy.constants import physical_constants as pc
from tqdm import tqdm

from vise.util.logger import get_logger
from vise.util.mix_in import ToCsvFileMixIn, ToJsonFileMixIn

logger = get_logger(__name__)

# Conversion factor from eV to cm^-1
EV_TO_INV_CM = pc["electron volt-inverse meter relationship"][0] / 100
eV_to_inv_cm = EV_TO_INV_CM  # Backward compatibility alias


def absorption_coeff(energy: float, real: float, imag: float) -> float:
    """Calculate absorption coefficient from dielectric function.

    Uses: α = 2√2π × √(√(ε₁² + ε₂²) - ε₁) × E × (conversion factor)

    Args:
        energy: Photon energy in eV.
        real: Real part of dielectric function (ε₁).
        imag: Imaginary part of dielectric function (ε₂).

    Returns:
        Absorption coefficient in cm⁻¹.
    """
    epsilon_magnitude = sqrt(real**2 + imag**2)
    return (
        2 * sqrt(2) * pi * sqrt(epsilon_magnitude - real)
        * energy * EV_TO_INV_CM
    )


def refractive_idx_real(e_real: float, e_imag: float) -> float:
    """Calculate real part of refractive index (n).

    n = √((ε₁ + √(ε₁² + ε₂²)) / 2)
    """
    epsilon_magnitude = sqrt(e_real**2 + e_imag**2)
    return sqrt(e_real + epsilon_magnitude) / sqrt(2)


def refractive_idx_imag(e_real: float, e_imag: float) -> float:
    """Calculate imaginary part of refractive index (extinction coefficient k).

    k = √((-ε₁ + √(ε₁² + ε₂²)) / 2)
    """
    epsilon_magnitude = sqrt(e_real**2 + e_imag**2)
    return sqrt(-e_real + epsilon_magnitude) / sqrt(2)


@dataclass
class DieleFuncData(MSONable, ToJsonFileMixIn, ToCsvFileMixIn):
    """Frequency-dependent dielectric function data.

    Stores real and imaginary parts of dielectric function along
    specified directions, with derived optical properties.

    Attributes:
        energies: Photon energies in eV.
        directions: Tensor directions (e.g., ["xx", "yy", "ave"]).
        diele_func_real: Real parts for each direction.
        diele_func_imag: Imaginary parts for each direction.
        band_gap: Band gap in eV (optional).
        title: Plot title (optional).
    """

    energies: List[float]
    directions: List[str]
    diele_func_real: List[List[float]]
    diele_func_imag: List[List[float]]
    band_gap: Optional[float] = None
    title: Optional[str] = None

    def __post_init__(self) -> None:
        """Validate data dimensions."""
        try:
            assert len(self.directions) == len(self.diele_func_real)
        except AssertionError:
            print(f"{len(self.directions)} vs {len(self.diele_func_real)}")
            raise

    @property
    def real_columns(self) -> List[str]:
        """Column names for real part data."""
        return [f"real_{d}" for d in self.directions]

    @property
    def imag_columns(self) -> List[str]:
        """Column names for imaginary part data."""
        return [f"imag_{d}" for d in self.directions]

    @property
    def absorption_columns(self) -> List[str]:
        """Column names for absorption coefficient data."""
        return [f"absorption_{d}" for d in self.directions]

    @property
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        data = {"energies(eV)": self.energies}

        for col_name, values in zip(self.real_columns, self.diele_func_real):
            data[col_name] = values
        for col_name, values in zip(self.imag_columns, self.diele_func_imag):
            data[col_name] = values
        for col_name, values in zip(self.absorption_columns, self.absorption_coeff):
            data[col_name] = values

        data["band_gap"] = [None] * len(self.energies)
        data["band_gap"][0] = self.band_gap

        return pd.DataFrame.from_dict(data)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "DieleFuncData":
        """Create from pandas DataFrame."""
        real_values: List[List[float]] = []
        imag_values: List[List[float]] = []
        directions: List[str] = []
        show_absorption_log = True

        # Handle both iteritems() (old pandas) and items() (new pandas)
        col_items = df.iteritems() if hasattr(df, "iteritems") else df.items()

        for column_name, item in col_items:
            if column_name in ["energies(eV)", "band_gap"]:
                continue
            elif "real" in column_name:
                real_values.append(item.tolist())
                directions.append(column_name.split("_")[-1])
            elif "imag" in column_name:
                imag_values.append(item.tolist())
            elif "absorption" in column_name:
                if show_absorption_log:
                    logger.info(
                        "Absorption column is skipped when constructing "
                        "DieleFuncData (derived from real/imag)."
                    )
                    show_absorption_log = False
            else:
                raise KeyError(f"Unknown column format: {column_name}")

        return cls(
            energies=df["energies(eV)"].tolist(),
            diele_func_real=real_values,
            diele_func_imag=imag_values,
            directions=directions,
            band_gap=float(df.loc[0, "band_gap"]),
        )

    @property
    def absorption_coeff(self) -> List[List[float]]:
        """Absorption coefficients in cm⁻¹ for each direction."""
        return [
            [
                absorption_coeff(e, r, i)
                for e, r, i in zip(self.energies, reals, imags)
            ]
            for reals, imags in zip(self.diele_func_real, self.diele_func_imag)
        ]

    @property
    def refractive_idx_real(self) -> List[List[float]]:
        """Real refractive index (n) for each direction."""
        return [
            [refractive_idx_real(r, i) for r, i in zip(reals, imags)]
            for reals, imags in zip(self.diele_func_real, self.diele_func_imag)
        ]

    @property
    def refractive_idx_imag(self) -> List[List[float]]:
        """Extinction coefficient (k) for each direction."""
        return [
            [refractive_idx_imag(r, i) for r, i in zip(reals, imags)]
            for reals, imags in zip(self.diele_func_real, self.diele_func_imag)
        ]

    @property
    def reflectivity(self) -> List[List[float]]:
        """Normal-incidence reflectivity R = ((n-1)² + k²) / ((n+1)² + k²)."""
        return [
            [
                ((n - 1) ** 2 + k ** 2) / ((n + 1) ** 2 + k ** 2)
                for n, k in zip(ns, ks)
            ]
            for ns, ks in zip(self.refractive_idx_real, self.refractive_idx_imag)
        ]


def min_e_w_target_coeff(
    energies: List[float],
    quantities: List[float],
    target_value: float,
) -> Optional[float]:
    """Find minimum energy where quantity exceeds target.

    Args:
        energies: Energy values.
        quantities: Corresponding quantity values.
        target_value: Target threshold.

    Returns:
        Energy at which quantity first exceeds target, or None.
    """
    for energy, quantity in zip(energies, quantities):
        if quantity > target_value:
            return energy

    logger.warning(
        f"Target value {target_value} not reached in the energy range."
    )
    return None


def make_shifted_diele_func(
    diele_func_data: DieleFuncData,
    original_band_gap: float,
    shift: float,
) -> DieleFuncData:
    """Create band-gap-shifted dielectric function using scissors operator.

    Args:
        diele_func_data: Original dielectric function data.
        original_band_gap: Original band gap in eV.
        shift: Energy shift to apply in eV.

    Returns:
        New DieleFuncData with shifted dielectric function.
    """
    shifted_imag = imag_shift(
        diele_func_data.diele_func_imag,
        diele_func_data.energies,
        original_band_gap + shift,
        shift,
    )
    shifted_real = kramers_kronig_trans(shifted_imag, diele_func_data.energies)

    return DieleFuncData(
        diele_func_data.energies,
        diele_func_data.directions,
        shifted_real.tolist(),
        shifted_imag.tolist(),
        original_band_gap + shift,
    )


def imag_shift(
    diele_func_imag: List[List[float]],
    energies: List[float],
    band_gap: float,
    shift: float,
) -> np.ndarray:
    """Apply scissors shift to imaginary dielectric function.

    Shifts the imaginary part to account for band gap correction,
    using linear interpolation.

    Args:
        diele_func_imag: Original imaginary dielectric function.
        energies: Energy grid.
        band_gap: New band gap after shift.
        shift: Energy shift amount.

    Returns:
        Shifted imaginary dielectric function.
    """
    energies_arr = np.array(energies)
    assert shift > 0

    result: List[List[float]] = []

    for energy_grid in energies_arr:
        old_energy = energy_grid - shift
        right_idx = np.argwhere(energies_arr > old_energy)[0][0]
        left_energy = energies_arr[right_idx - 1]
        right_energy = energies_arr[right_idx]

        # Linear interpolation weight
        left_weight = (right_energy - old_energy) / (right_energy - left_energy)

        inner_result: List[float] = []
        for imag_idx in range(len(diele_func_imag)):
            if energy_grid < band_gap:
                inner_result.append(0.0)
            else:
                # Interpolate from original energies
                interpolated = (
                    diele_func_imag[imag_idx][right_idx - 1] * left_weight
                    + diele_func_imag[imag_idx][right_idx] * (1 - left_weight)
                )
                # Scale by energy ratio
                scaled = interpolated * (energy_grid - shift) / energy_grid
                inner_result.append(scaled)

        result.append(inner_result)

    return np.array(result).T


def kramers_kronig_trans(
    diele_func_imag: np.ndarray,
    energies: List[float],
    ita: float = 0.01,
) -> np.ndarray:
    """Apply Kramers-Kronig transformation to get real part.

    Computes real dielectric function from imaginary part using
    the Kramers-Kronig relation.

    Args:
        diele_func_imag: Imaginary dielectric function [direction, energy].
        energies: Energy grid.
        ita: Small broadening parameter to avoid singularities.

    Returns:
        Real dielectric function.
    """
    mesh = energies[1] - energies[0]
    result: List[List[float]] = []

    # Precompute energy differences
    ee2ss = [
        [e**2 - energy_grid**2 for e in energies]
        for energy_grid in energies
    ]

    for imag_idx in tqdm(range(len(diele_func_imag))):
        imags = diele_func_imag[imag_idx, :]

        # Reuse previous result if identical
        if imag_idx > 0 and np.allclose(imags, diele_func_imag[imag_idx - 1, :]):
            pass  # inner_result unchanged from previous iteration
        elif np.count_nonzero(imags) == 0:
            inner_result = [0.0] * len(energies)
        else:
            inner_result = []
            for ee2s in ee2ss:
                integrals = [
                    e * imag * ee2 / (ee2**2 + ita**2)
                    for e, ee2, imag in zip(energies, ee2s, imags)
                ]
                integral = sum(integrals) * mesh * 2 / pi

                # Add 1 for diagonal components (ε_∞ = 1)
                if imag_idx < 3:
                    integral += 1

                inner_result.append(integral)

        result.append(inner_result)

    return np.array(result)
