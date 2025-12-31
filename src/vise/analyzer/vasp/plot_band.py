# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Band structure plot data extraction from VASP.

This module provides utilities for creating band plot information
from VASP band structure calculations.
"""

import re
from copy import deepcopy
from typing import Dict, List, Optional

import numpy as np
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp import Vasprun

from vise.analyzer.plot_band import (
    BandEdgeForPlot,
    BandEnergyInfo,
    BandPlotInfo,
    XTicks,
)
from vise.analyzer.plot_brillouin_zone import BZPlotInfo
from vise.util.string import latexify


def greek_to_unicode(label: str) -> str:
    """Convert Greek letter names to Unicode symbols.

    Args:
        label: K-point label (e.g., "GAMMA", "GM").

    Returns:
        Label with Greek characters (e.g., "Γ").

    Note:
        Order matters as SIGMA contains GM substring.
    """
    replacements = {
        "GAMMA": "Γ",
        "SIGMA": "Σ",
        "GM": "Γ",
        "DELTA": "Δ",
        "LAMBDA": "Λ",
    }
    for old, new in replacements.items():
        label = label.replace(old, new)
    return label


def italic_to_roman(label: str) -> str:
    """Convert italic subscripted labels to roman font in LaTeX.

    Args:
        label: Label with subscripts (e.g., "S_0").

    Returns:
        Label with roman font subscripts (e.g., "{\\rm S}_0").
    """
    return re.sub(r"([A-Z])_([0-9])", r"{\\rm \1}_\2", label)


def _sanitize_label(label: str) -> str:
    """Apply all label sanitizations."""
    return italic_to_roman(greek_to_unicode(label))


def _sanitize_labels(labels: List[str]) -> List[str]:
    """Apply sanitization to a list of labels."""
    return [_sanitize_label(label) for label in labels]


class BandPlotInfoFromVasp:
    """Extract band plot information from VASP calculations.

    Supports comparison of two band structures (e.g., PBE vs GW)
    with the same k-path.

    Attributes:
        bs: Primary band structure.
        bs2: Optional secondary band structure for comparison.
    """

    def __init__(
        self,
        vasprun: Vasprun,
        kpoints_filename: str,
        first_band_plot_name: Optional[str] = None,
        wavecar_name: Optional[str] = None,
        vasprun2: Optional[Vasprun] = None,
        second_band_plot_name: Optional[str] = None,
        energy_window: Optional[List[float]] = None,
    ) -> None:
        """Initialize from VASP output files.

        Args:
            vasprun: Primary vasprun.xml.
            kpoints_filename: Path to KPOINTS file.
            first_band_plot_name: Label for first band structure.
            wavecar_name: Path to WAVECAR (for irreps, currently unused).
            vasprun2: Optional second vasprun.xml for comparison.
            second_band_plot_name: Label for second band structure.
            energy_window: [min, max] energy range to include bands.
        """
        self._composition = vasprun.final_structure.composition
        self.bs = vasprun.get_band_structure(kpoints_filename, line_mode=True)
        self.first_band_plot_name = first_band_plot_name or "1"
        self.wavecar_name = wavecar_name
        self.rlat = vasprun.final_structure.lattice.reciprocal_lattice
        self.kpoints = self.bs.kpoints
        self.energy_window = energy_window

        # Handle optional second band structure
        if vasprun2:
            assert vasprun.final_structure == vasprun2.final_structure
            self.vasprun2 = vasprun2
            self.bs2 = vasprun2.get_band_structure(kpoints_filename, line_mode=True)
            self.second_band_plot_name = second_band_plot_name or "2"
        else:
            self.vasprun2 = None
            self.bs2 = None
            self.second_band_plot_name = None

    def make_band_plot_info(self) -> BandPlotInfo:
        """Create BandPlotInfo from VASP data.

        Returns:
            BandPlotInfo ready for plotting.
        """
        bs_plotter = BSPlotter(self.bs)
        plot_data = bs_plotter.bs_plot_data(zero_to_efermi=False)
        distances = [list(d) for d in plot_data["distances"]]

        # Create primary band info
        band_info_1 = BandEnergyInfo(
            band_energies=self._remove_spin_key(plot_data),
            band_edge=self._band_edge(self.bs, plot_data),
            fermi_level=self.bs.efermi,
        )
        band_infos = {self.first_band_plot_name: band_info_1}

        # Add secondary band structure if provided
        if self.vasprun2:
            plot_data2 = BSPlotter(self.bs2).bs_plot_data(zero_to_efermi=False)
            band_infos[self.second_band_plot_name] = BandEnergyInfo(
                band_energies=self._remove_spin_key(plot_data2),
                band_edge=self._band_edge(self.bs2, plot_data2),
                fermi_level=self.bs2.efermi,
            )

        # Extract x-axis tick information
        ticks = bs_plotter.get_ticks_old()
        x_ticks = XTicks(_sanitize_labels(ticks["label"]), ticks["distance"])

        return BandPlotInfo(
            band_energy_infos=band_infos,
            distances_by_branch=distances,
            x_ticks=x_ticks,
            title=self._title,
        )

    def make_bz_plot_info(self) -> BZPlotInfo:
        """Create Brillouin zone plot information.

        Returns:
            BZPlotInfo with zone faces, labels, and band paths.
        """
        # Get Wigner-Seitz cell faces
        faces = [
            [[float(k) for k in j] for j in i]
            for i in self.rlat.get_wigner_seitz_cell()
        ]

        labels: Dict[str, Dict[str, List[float]]] = {}
        band_paths: List[List[List[float]]] = []
        concat = False
        init_point: Optional[List[float]] = None

        for kpoint in self.kpoints:
            if kpoint.label:
                cart_coords = list(kpoint.cart_coords)
                frac_coords = list(kpoint.frac_coords)
                label = greek_to_unicode(kpoint.label)
                labels[label] = {"cart": cart_coords, "frac": frac_coords}

                # Track band paths between labeled points
                if not concat and init_point:
                    band_paths.append([init_point, cart_coords])
                init_point = cart_coords
                concat = True
            else:
                concat = False

        return BZPlotInfo(faces, labels, band_paths, self.rlat.matrix.tolist())

    def _remove_spin_key(self, plot_data: dict) -> List[List[List[List[float]]]]:
        """Reshape energy data from pymatgen format.

        Pymatgen format: {Spin: [np.array(n_bands, n_kpts), ...]}
        Output format: [branch][spin][band][kpoint]
        """
        num_spin = len(plot_data["energy"])
        num_branch = len(plot_data["energy"]["1"])

        result = [[[] for _ in range(num_spin)] for _ in range(num_branch)]

        # Sort by spin key (descending so spin up comes first)
        for spin_idx, (_, branch_energies) in enumerate(
            sorted(plot_data["energy"].items(), key=lambda x: x[0], reverse=True)
        ):
            for branch_idx, branch_energy in enumerate(branch_energies):
                if self.energy_window:
                    # Filter bands by energy window
                    keep_indices = []
                    for i in range(len(branch_energy)):
                        band_max = np.max(branch_energy[i, :])
                        band_min = np.min(branch_energy[i, :])
                        if self._in_energy_window(band_max, band_min):
                            keep_indices.append(i)
                    filtered = np.delete(
                        branch_energy,
                        [i for i in range(len(branch_energy)) if i not in keep_indices],
                        axis=0,
                    )
                    energies = filtered.tolist()
                else:
                    energies = branch_energy.tolist()

                result[branch_idx][spin_idx] = deepcopy(energies)

        return result

    def _in_energy_window(self, band_max: float, band_min: float) -> bool:
        """Check if band overlaps with energy window."""
        return (
            band_max >= self.energy_window[0] and band_min <= self.energy_window[1]
        )

    @staticmethod
    def _band_edge(
        bs: BandStructure, plot_data: dict
    ) -> Optional[BandEdgeForPlot]:
        """Extract band edge information for semiconductors."""
        if bs.is_metal():
            return None

        return BandEdgeForPlot(
            vbm=plot_data["vbm"][0][1],
            cbm=plot_data["cbm"][0][1],
            vbm_distances=[i[0] for i in plot_data["vbm"]],
            cbm_distances=[i[0] for i in plot_data["cbm"]],
        )

    @property
    def _title(self) -> str:
        """Generate plot title from composition."""
        return latexify(self._composition.reduced_formula)
