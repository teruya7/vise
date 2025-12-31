# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""VASP-specific DOS data extraction.

This module provides utilities for extracting DOS data from VASP
vasprun.xml output files.
"""

from typing import List, Optional

import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Vasprun

from vise.analyzer.dos_data import DosData, PDos
from vise.util.logger import get_logger

logger = get_logger(__name__)


class DosDataFromVasp:
    """DOS data extractor for VASP calculations.

    Extracts total and partial DOS from vasprun.xml, optionally
    cropping to a specified energy window.

    Attributes:
        complete_dos: Complete DOS from pymatgen.
        vertical_lines: Energies for vertical markers.
        base_energy: Reference energy for alignment.
        energy_window: Optional energy range filter.
    """

    def __init__(
        self,
        vasprun: Vasprun,
        vertical_lines: Optional[List[float]] = None,
        base_energy: float = 0.0,
        crop_first_value: bool = False,
        energy_window: Optional[List[float]] = None,
    ) -> None:
        """Initialize the DOS extractor.

        Args:
            vasprun: Parsed vasprun.xml.
            vertical_lines: Energies for vertical markers (VBM, CBM, etc.).
            base_energy: Reference energy for alignment.
            crop_first_value: Whether to skip the first energy point
                            (sometimes contains noise).
            energy_window: [min, max] energy range to include.
        """
        self.complete_dos = vasprun.complete_dos
        self.vertical_lines = vertical_lines or []
        self.base_energy = base_energy
        self.energy_window = energy_window

        self.min_energy_idx = 1 if crop_first_value else 0
        self.max_energy_idx = len(self.complete_dos.energies)

    def make_dos_data(self) -> DosData:
        """Create DosData from VASP output.

        Returns:
            DosData object with total and partial DOS.
        """
        energies = self.complete_dos.energies.tolist()

        # Apply energy window filter if specified
        if self.energy_window:
            min_idx = self._get_min_energy_idx(energies)
            max_idx = self._get_max_energy_idx(energies)
            if min_idx:
                self.min_energy_idx = min_idx
            if max_idx:
                self.max_energy_idx = max_idx

        energies = energies[self.min_energy_idx : self.max_energy_idx + 1]

        return DosData(
            energies=energies,
            total=np.array(self._total),
            pdos=self._pdos,
            vertical_lines=self.vertical_lines,
            base_energy=self.base_energy,
        )

    def _get_max_energy_idx(self, energies: List[float]) -> Optional[int]:
        """Find index of maximum energy in window."""
        for i, energy in enumerate(energies):
            if energy > self.energy_window[1]:
                return i - 1
        return None

    def _get_min_energy_idx(self, energies: List[float]) -> Optional[int]:
        """Find index of minimum energy in window."""
        for i, energy in enumerate(energies):
            if energy >= self.energy_window[0]:
                return i
        return None

    @property
    def _pdos(self) -> List[PDos]:
        """Extract partial DOS for each atom."""
        result: List[PDos] = []

        for site_dos in self.complete_dos.pdos.values():
            pdos_kwargs = {}

            for orbital, orbital_dos in site_dos.items():
                # Extract DOS for each spin channel
                pdos = np.array([
                    orbital_dos[spin]
                    for spin in [Spin.up, Spin.down]
                    if spin in orbital_dos
                ])

                # Apply energy window slice
                pdos_kwargs[str(orbital)] = pdos[
                    :, self.min_energy_idx : self.max_energy_idx + 1
                ]

            try:
                result.append(PDos.from_dict(pdos_kwargs))
            except TypeError:
                logger.warning(
                    "Orbital-decomposed DOS required. Set LORBIT = 11 in VASP."
                )
                raise

        return result

    @property
    def _total(self) -> List[np.ndarray]:
        """Extract total DOS for each spin channel."""
        result: List[np.ndarray] = []

        for spin in [Spin.up, Spin.down]:
            if spin in self.complete_dos.densities:
                dos = self.complete_dos.densities[spin]
                result.append(dos[self.min_energy_idx : self.max_energy_idx + 1])

        return result
