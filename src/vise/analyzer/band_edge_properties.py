# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Band edge property analysis.

This module provides classes and functions for analyzing valence band
maximum (VBM) and conduction band minimum (CBM) properties from
electronic structure calculations.
"""

from copy import deepcopy
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from monty.json import MSONable
from numpy import argmin
from pymatgen.electronic_structure.core import Spin

from vise.defaults import defaults
from vise.util.logger import get_logger

logger = get_logger(__name__)


@dataclass
class BandEdge(MSONable):
    """Information about a band edge (VBM or CBM).

    Attributes:
        energy: Energy of the band edge in eV.
        spin: Spin channel (up or down).
        band_index: Index of the band at this edge.
        kpoint_coords: K-point coordinates in reciprocal fractional coords.
        kpoint_index: Index of the k-point in the mesh.
        data_source: Source of data (e.g., 'vasp', 'materials_project').
        symbol: High-symmetry point symbol (e.g., 'Γ', 'X').
    """

    energy: float
    spin: Spin
    band_index: int
    kpoint_coords: List[float]
    kpoint_index: Optional[int] = None
    data_source: Optional[str] = None
    symbol: Optional[str] = None

    def is_direct(self, other: "BandEdge") -> bool:
        """Check if this edge forms a direct gap with another.

        Args:
            other: The other band edge (CBM if this is VBM, or vice versa).

        Returns:
            True if both edges occur at the same spin and k-point.
        """
        return (
            self.spin == other.spin and
            self.kpoint_coords == other.kpoint_coords
        )

    def __repr__(self) -> str:
        kpt_str = (
            f"{self.kpoint_coords[0]:5.3f} "
            f"{self.kpoint_coords[1]:5.3f} "
            f"{self.kpoint_coords[2]:5.3f}"
        )
        return (
            f"energy position: {self.energy}, spin: {self.spin.name:>4}, "
            f"band index {self.band_index}, "
            f"k-point index {self.kpoint_index}, k-point coords {kpt_str}"
        )

    def as_dict(self) -> dict:
        """Serialize to dictionary for JSON."""
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "energy": self.energy,
            "spin": int(self.spin),
            "band_index": self.band_index,
            "kpoint_index": self.kpoint_index,
            "kpoint_coords": self.kpoint_coords,
            "data_source": self.data_source,
            "symbol": self.symbol,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "BandEdge":
        """Create from dictionary."""
        kwargs = deepcopy(d)
        kwargs["spin"] = Spin(d["spin"])
        kwargs.pop("@module", None)
        kwargs.pop("@class", None)
        return cls(**kwargs)


class BandEdgeProperties:
    """Analyzer for band edge properties from eigenvalue data.

    Determines VBM and CBM positions, band gap, and whether the gap
    is direct or indirect.

    Attributes:
        vbm_info: BandEdge information for VBM (None if metal).
        cbm_info: BandEdge information for CBM (None if metal).
    """

    def __init__(
        self,
        eigenvalues: Dict[Spin, np.ndarray],
        nelect: float,
        magnetization: float,
        kpoint_coords: List[List[float]],
        integer_criterion: float = defaults.integer_criterion,
        is_non_collinear: bool = False,
    ) -> None:
        """Initialize with eigenvalue data.

        Args:
            eigenvalues: Eigenvalues indexed [spin][kpoint][band].
            nelect: Total number of electrons.
            magnetization: Total magnetization in Bohr magnetons.
            kpoint_coords: K-point coordinates [[k1_xyz], [k2_xyz], ...].
            integer_criterion: Threshold for determining integer occupation.
            is_non_collinear: Whether calculation is non-collinear.
        """
        assert 0 < integer_criterion < 0.5

        self._eigenvalues = eigenvalues
        self._nelect = nelect
        self._magnetization = magnetization  # In Bohr magneton
        self._kpoint_coords = kpoint_coords
        self._integer_criterion = integer_criterion
        self._is_non_collinear = is_non_collinear

        self._calculate_vbm_cbm()

        if self._is_metal():
            self.vbm_info = None
            self.cbm_info = None

    def _calculate_vbm_cbm(self) -> None:
        """Find VBM and CBM from eigenvalues."""
        vbm_energy = float("-inf")
        cbm_energy = float("inf")

        for spin, eigenvalues in self._eigenvalues.items():
            ho_band_idx = self._ho_band_index(spin)
            lu_band_idx = ho_band_idx + 1

            # Find highest occupied (VBM candidate)
            ho_eigenvalue = np.amax(eigenvalues[:, ho_band_idx])
            if ho_eigenvalue > vbm_energy:
                self.vbm_info = self.band_edge(
                    eigenvalues, ho_band_idx, ho_eigenvalue, spin
                )
                vbm_energy = ho_eigenvalue

            # Find lowest unoccupied (CBM candidate)
            lu_eigenvalue = np.amin(eigenvalues[:, lu_band_idx])
            if lu_eigenvalue < cbm_energy:
                self.cbm_info = self.band_edge(
                    eigenvalues, lu_band_idx, lu_eigenvalue, spin
                )
                cbm_energy = lu_eigenvalue

    def band_edge(
        self,
        eigenvalues: np.ndarray,
        band_index: int,
        eigenvalue: float,
        spin: Spin,
    ) -> BandEdge:
        """Create BandEdge from eigenvalue data."""
        k_index = int(np.where(eigenvalues[:, band_index] == eigenvalue)[0][0])
        return BandEdge(
            float(eigenvalue),
            spin,
            band_index,
            self._kpoint_coords[k_index],
            k_index,
        )

    def _ho_band_index(self, spin: Spin) -> int:
        """Get highest occupied band index for given spin."""
        if self._is_non_collinear:
            return round(self._nelect) - 1

        if spin == Spin.up:
            num_occupied = (self._nelect + self._magnetization) / 2
        else:
            num_occupied = (self._nelect - self._magnetization) / 2

        return round(num_occupied) - 1

    def _is_metal(self) -> bool:
        """Check if system is metallic."""
        nelect_frac = abs(round(self._nelect) - self._nelect)
        is_nelect_frac = nelect_frac > self._integer_criterion

        mag_frac = abs(round(self._magnetization) - self._magnetization)
        is_mag_frac = mag_frac > self._integer_criterion

        is_vbm_higher = self.vbm_info.energy > self.cbm_info.energy

        return is_nelect_frac or is_mag_frac or is_vbm_higher

    @property
    def is_metal(self) -> bool:
        """Whether the system is metallic."""
        return self.vbm_info is None

    @property
    def is_direct(self) -> Optional[bool]:
        """Whether the band gap is direct."""
        if self.vbm_info:
            return self.vbm_info.is_direct(self.cbm_info)
        return None

    @property
    def band_gap(self) -> Optional[float]:
        """Band gap in eV (None for metals)."""
        if self.vbm_info:
            return self.cbm_info.energy - self.vbm_info.energy
        return None

    def min_gap_w_coords(
        self, same_energy_threshold: float = 1e-5
    ) -> Tuple[float, List[List[float]]]:
        """Find minimum direct gap and its k-point positions.

        Note that there may be multiple k-points with the same minimum gap.

        Args:
            same_energy_threshold: Energy tolerance for identifying
                equivalent gap positions.

        Returns:
            Tuple of (min_gap, list of k-point coordinates).
        """
        min_gap = float("inf")
        kpoint_indices: List[int] = []

        for spin, eigenvalues in self._eigenvalues.items():
            ho_idx = self._ho_band_index(spin)
            lu_idx = ho_idx + 1

            ho_eigs = eigenvalues[:, ho_idx]
            lu_eigs = eigenvalues[:, lu_idx]
            direct_gaps = lu_eigs - ho_eigs

            min_k = argmin(direct_gaps)
            if direct_gaps[min_k] < min_gap + same_energy_threshold:
                min_gap = direct_gaps[min_k]
                kpoint_indices.append(min_k)

        kpoint_coords = [self._kpoint_coords[i] for i in kpoint_indices]
        return min_gap, kpoint_coords

    @property
    def vbm_cbm(self) -> Optional[List[float]]:
        """VBM and CBM energies as [vbm, cbm]."""
        if self.vbm_info:
            return [self.vbm_info.energy, self.cbm_info.energy]
        return None

    def __repr__(self) -> str:
        if self.vbm_info:
            lines = [
                f"Band gap {self.band_gap:5.3f} eV",
                f"VBM {self.vbm_info}",
                f"CBM {self.cbm_info}",
            ]
            return "\n".join(lines)
        return "Metal"


def is_band_gap(
    band_gap: Optional[float],
    vbm_cbm: Optional[List[float]],
    show_info: bool = False,
) -> bool:
    """Check if system has a band gap above threshold.

    Args:
        band_gap: Direct band gap value (if known).
        vbm_cbm: [VBM, CBM] energies (used if band_gap not provided).
        show_info: Whether to log the result.

    Returns:
        True if band gap exists and exceeds threshold.
    """
    if not band_gap:
        band_gap = vbm_cbm[1] - vbm_cbm[0] if vbm_cbm else None

    if band_gap:
        has_gap = band_gap > defaults.band_gap_criterion
        if show_info:
            comparison = ">" if has_gap else "<"
            logger.info(
                f"Is there a band gap?: {has_gap} "
                f"({round(band_gap, 3)} {comparison} "
                f"{defaults.band_gap_criterion} eV)"
            )
        return has_gap

    return False


def merge_band_edge(
    be_1: BandEdge,
    be_2: BandEdge,
    edge: str,
    threshold: float = 0.001,
) -> BandEdge:
    """Merge two BandEdge objects, keeping the more extreme one.

    Args:
        be_1: First band edge.
        be_2: Second band edge.
        edge: 'vbm' (keep higher energy) or 'cbm' (keep lower energy).
        threshold: Energy tolerance for considering edges equivalent.

    Returns:
        Merged BandEdge with combined data_source if equivalent.

    Raises:
        ValueError: If edge is not 'vbm' or 'cbm'.
    """
    if edge not in ["vbm", "cbm"]:
        raise ValueError("edge must be 'vbm' or 'cbm'")

    if edge == "vbm":
        result = be_1 if be_1.energy > be_2.energy else be_2
    else:
        result = be_1 if be_1.energy < be_2.energy else be_2

    # Merge data sources if edges are equivalent
    if (
        abs(be_1.energy - be_2.energy) < threshold
        and np.allclose(be_1.kpoint_coords, be_2.kpoint_coords)
        and be_1.spin == be_2.spin
    ):
        result.data_source = f"{be_1.data_source} {be_2.data_source}"

    return deepcopy(result)
