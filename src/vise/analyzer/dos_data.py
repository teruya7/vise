# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Density of states data structures.

This module provides data classes for representing and manipulating
density of states (DOS) data, including partial DOS (PDOS) decomposed
by orbital character.
"""

from collections import defaultdict
from copy import copy, deepcopy
from dataclasses import dataclass
from functools import reduce
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from monty.json import MSONable

from vise.util.logger import get_logger
from vise.util.mix_in import ToCsvFileMixIn, ToJsonFileMixIn

logger = get_logger(__name__)


@dataclass
class PDos(MSONable, ToJsonFileMixIn):
    """Partial density of states decomposed by orbital.

    Stores DOS contributions from each angular momentum channel
    (s, p, d, f) and their m_l components.

    Attributes:
        s: s-orbital DOS [spin][energy].
        px, py, pz: p-orbital components.
        dxy, dyz, dxz, dx2, dz2: d-orbital components.
        f_3 through f3: f-orbital components (optional).
    """

    s: np.ndarray
    px: np.ndarray
    py: np.ndarray
    pz: np.ndarray
    dxy: np.ndarray
    dyz: np.ndarray
    dxz: np.ndarray
    dx2: np.ndarray
    dz2: np.ndarray
    f_3: Optional[np.ndarray] = None
    f_2: Optional[np.ndarray] = None
    f_1: Optional[np.ndarray] = None
    f0: Optional[np.ndarray] = None
    f1: Optional[np.ndarray] = None
    f2: Optional[np.ndarray] = None
    f3: Optional[np.ndarray] = None

    @property
    def p(self) -> np.ndarray:
        """Total p-orbital DOS."""
        return self.px + self.py + self.pz

    @property
    def d(self) -> np.ndarray:
        """Total d-orbital DOS."""
        return self.dxy + self.dyz + self.dxz + self.dx2 + self.dz2

    @property
    def f(self) -> Optional[np.ndarray]:
        """Total f-orbital DOS (None if not available)."""
        try:
            return (
                self.f_3 + self.f_2 + self.f_1
                + self.f0 + self.f1 + self.f2 + self.f3
            )
        except TypeError:
            return None

    @classmethod
    def from_dict(cls, d: dict) -> "PDos":
        """Create PDos from dictionary.

        Handles both fully decomposed (px, py, pz) and partially
        decomposed (p, d, f) input formats.
        """
        if "px" in d:
            return super().from_dict(d)

        # Handle non-decomposed format by equal division
        logger.warning(
            "PDOS not decomposed to individual orbitals (px, py, pz). "
            "Dividing equally among components."
        )
        result = {
            "s": d["s"],
            "px": d["p"] / 3,
            "py": d["p"] / 3,
            "pz": d["p"] / 3,
            "dxy": d["d"] / 5,
            "dyz": d["d"] / 5,
            "dxz": d["d"] / 5,
            "dx2": d["d"] / 5,
            "dz2": d["d"] / 5,
        }
        if "f" in d:
            result.update({
                "f_3": d["f"] / 7,
                "f_2": d["f"] / 7,
                "f_1": d["f"] / 7,
                "f0": d["f"] / 7,
                "f1": d["f"] / 7,
                "f2": d["f"] / 7,
                "f3": d["f"] / 7,
            })
        return cls(**result)

    def __add__(self, other: "PDos") -> "PDos":
        """Sum two PDos objects element-wise."""
        args = {}
        for key, value in self.__dict__.items():
            if value is None:
                args[key] = None
            else:
                args[key] = value + getattr(other, key)
        return PDos(**args)

    def __eq__(self, other: "PDos") -> bool:
        """Check equality with another PDos."""
        for key, value in self.__dict__.items():
            if (value != getattr(other, key)).any():
                return False
        return True


@dataclass
class DosData(MSONable):
    """Complete DOS data including total and partial contributions.

    Attributes:
        energies: Energy grid in eV.
        total: Total DOS [spin][energy].
        pdos: Partial DOS for each atom.
        vertical_lines: Energies for vertical markers (e.g., VBM, CBM).
        base_energy: Reference energy (e.g., Fermi level or VBM).
    """

    energies: List[float]
    total: np.ndarray
    pdos: List[PDos]
    vertical_lines: List[float]
    base_energy: Optional[float] = 0.0

    def __post_init__(self) -> None:
        """Ensure total is numpy array."""
        if isinstance(self.total, list):
            self.total = np.array(self.total)

    @property
    def spin(self) -> bool:
        """Whether data is spin-polarized."""
        return len(self.total) != 1

    def dos_plot_data(
        self,
        grouped_atom_indices: Dict[str, List[int]],
        energy_range: Optional[List[float]] = None,
        dos_ranges: Optional[List[List[float]]] = None,
        title: Optional[str] = None,
    ) -> "DosPlotData":
        """Create DosPlotData for plotting.

        Args:
            grouped_atom_indices: Atom groups keyed by name.
            energy_range: Energy range for plot [min, max].
            dos_ranges: Y-axis ranges for each panel.
            title: Plot title.

        Returns:
            DosPlotData ready for plotting.
        """
        if dos_ranges is not None:
            print(len(grouped_atom_indices))
            print(len(dos_ranges))
            try:
                # total + pdos panels
                assert len(grouped_atom_indices) + 1 == len(dos_ranges)
            except AssertionError:
                print(
                    f"Number of DOS ranges ({len(dos_ranges)}) doesn't match "
                    f"number of panels ({len(grouped_atom_indices) + 1}). "
                    f"Note: total DOS is included as first panel."
                )
                raise

        # Build doses list starting with total
        doses: List[List[DosBySpinEnergy]] = [
            [DosBySpinEnergy("", self.total.tolist())]
        ]
        names = ["total"]

        for name, atom_indices in grouped_atom_indices.items():
            # Sum PDos for all atoms in group
            pdos_list = [self.pdos[idx] for idx in atom_indices]
            pdos = reduce(lambda x, y: x + y, pdos_list)

            orbital_doses = [
                DosBySpinEnergy("s", pdos.s.tolist()),
                DosBySpinEnergy("p", pdos.p.tolist()),
                DosBySpinEnergy("d", pdos.d.tolist()),
            ]
            if pdos.f is not None:
                orbital_doses.append(DosBySpinEnergy("f", pdos.f.tolist()))

            doses.append(orbital_doses)
            names.append(name)

        energy_range = energy_range or [-5, 10]
        abs_xlim = [x + self.base_energy for x in energy_range]

        if dos_ranges is None:
            dos_ranges = default_dos_ranges(
                abs_xlim, doses, self.energies, self.spin
            )

        # Shift energies relative to base
        relative_energies = [e - self.base_energy for e in self.energies]
        shifted_lines = [e - self.base_energy for e in self.vertical_lines]

        return DosPlotData(
            relative_energies, doses, names, energy_range,
            dos_ranges, shifted_lines, title
        )


def default_dos_ranges(
    energy_range: List[float],
    doses: List[List["DosBySpinEnergy"]],
    energies: List[float],
    spin_polarized: bool,
    multi: float = 1.1,
    round_digit: int = 2,
) -> List[List[float]]:
    """Calculate default y-axis ranges for DOS plots.

    Args:
        energy_range: Energy range [min, max].
        doses: DOS data for each panel.
        energies: Energy grid.
        spin_polarized: Whether data is spin-polarized.
        multi: Multiplier for max DOS to add padding.
        round_digit: Decimal places for rounding.

    Returns:
        List of [y_min, y_max] for each panel.
    """
    mask = np.ma.masked_outside(energies, energy_range[0], energy_range[1]).mask

    max_dos_by_panel = []
    for panel_dos in doses:
        max_dos_by_panel.append(
            np.max([dos.max_dos(mask) for dos in panel_dos])
        )

    total_dos_max = max_dos_by_panel[0]
    plot_maxes = [total_dos_max]

    if len(max_dos_by_panel) > 1:
        pdos_max = max(max_dos_by_panel[1:])
        plot_maxes.extend([pdos_max] * (len(doses) - 1))

    max_y_ranges = [round(val * multi, round_digit) for val in plot_maxes]

    if spin_polarized:
        return [[-y, y] for y in max_y_ranges]
    return [[0, y] for y in max_y_ranges]


@dataclass
class DosBySpinEnergy(MSONable):
    """DOS data for a single orbital type.

    Attributes:
        name: Orbital name (e.g., "s", "p", "d", "" for total).
        dos: DOS values [spin][energy].
    """

    name: str
    dos: List[List[float]]

    def max_dos(self, mask: Optional[List[bool]] = None) -> float:
        """Get maximum DOS value within mask."""
        return max(
            np.max(np.ma.masked_array(d, mask).compressed())
            for d in self.dos
        )


@dataclass
class DosPlotData(MSONable, ToJsonFileMixIn, ToCsvFileMixIn):
    """Data prepared for DOS plotting.

    Attributes:
        relative_energies: Energies relative to reference.
        doses: DOS for each panel [panel][orbital].
        names: Name for each panel.
        energy_range: X-axis range.
        dos_ranges: Y-axis range for each panel.
        energy_lines: Vertical line positions (VBM, CBM, etc.).
        title: Plot title.
    """

    relative_energies: List[float]
    doses: List[List[DosBySpinEnergy]]
    names: List[str]
    energy_range: List[float]
    dos_ranges: List[List[float]]
    energy_lines: List[float]
    title: Optional[str] = None

    @property
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        columns = ["energy(eV)"]
        inv_data = [list(self.relative_energies)]

        for name, panel_dos in zip(self.names, self.doses):
            for dos_by_orbital in panel_dos:
                for spin, spin_dos in zip(["up", "down"], dos_by_orbital.dos):
                    columns.append(f"{name}_{dos_by_orbital.name}_{spin}")
                    inv_data.append(spin_dos)

        columns.append("energy_lines")
        inv_data.append(self.energy_lines)

        df = pd.DataFrame(inv_data).T
        df.columns = columns
        return df

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "DosPlotData":
        """Create from pandas DataFrame."""
        sanitized_data = defaultdict(lambda: defaultdict(list))

        for full_name in df.columns[1:-1]:
            name, orbital, spin = full_name.split("_")
            sanitized_data[name][orbital].append(list(df[full_name]))

        doses, names = [], []
        for key, orbitals in sanitized_data.items():
            names.append(key)
            doses.append([
                {"name": orbital_name, "dos": orbital_dos}
                for orbital_name, orbital_dos in orbitals.items()
            ])

        d = {
            "relative_energies": df["energy(eV)"].tolist(),
            "doses": doses,
            "names": names,
            "energy_lines": [i for i in df.energy_lines if i],
        }
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls, d: dict) -> "DosPlotData":
        """Create from dictionary with backward compatibility."""
        for key in copy(d):
            if key.startswith("@"):
                d.pop(key)

        # Backward compatibility for old field names
        if "xlim" in d:
            d["energy_range"] = d.pop("xlim")
        if "ylim_set" in d:
            d["dos_ranges"] = d.pop("ylim_set")
        if "vertical_lines" in d:
            d["energy_lines"] = d.pop("vertical_lines")

        # Convert nested dicts to DosBySpinEnergy
        for i, panel in enumerate(d["doses"]):
            for j, orbital in enumerate(panel):
                d["doses"][i][j] = DosBySpinEnergy.from_dict(orbital)

        # Set defaults if missing
        if "energy_range" not in d:
            d["energy_range"] = [-5, 10]
        if "dos_ranges" not in d:
            spin_polarized = len(d["doses"][0][0].dos) == 2
            d["dos_ranges"] = default_dos_ranges(
                d["energy_range"],
                d["doses"],
                d["relative_energies"],
                spin_polarized,
            )

        return cls(**d)


def scissor_energy(
    dos_plot_data: DosPlotData,
    energy_shift: float,
) -> DosPlotData:
    """Apply scissors operator to shift conduction band energies.

    Shifts energies above the band gap middle by the specified amount.

    Args:
        dos_plot_data: Original DOS plot data.
        energy_shift: Energy to add to conduction band states.

    Returns:
        New DosPlotData with shifted energies.
    """
    result = deepcopy(dos_plot_data)
    gap_middle = np.mean(dos_plot_data.energy_lines)

    # Find where energies exceed gap middle
    above_gap = np.array(result.relative_energies) >= gap_middle
    shift_start_idx = int(np.argmax(above_gap))

    # Shift energies above gap
    for i in range(shift_start_idx, len(result.relative_energies)):
        result.relative_energies[i] += energy_shift

    # Shift CBM energy line
    result.energy_lines[1] += energy_shift

    return result
