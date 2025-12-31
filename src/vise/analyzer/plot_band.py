# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Band structure plotting utilities.

This module provides data classes and plotters for electronic band
structure visualization, including irreducible representations and
support for multiple band structure comparisons.
"""

from copy import deepcopy
from dataclasses import dataclass
from itertools import cycle
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from monty.json import MSONable

from vise.error import ViseError
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

logger = get_logger(__name__)


@dataclass
class Irrep(MSONable):
    """Irreducible representation at a k-point.

    Stores symmetry classification of bands at a high-symmetry point.

    Attributes:
        frac_coords: K-point coordinates in reciprocal fractional coords.
        symbols: Irrep symbols for each band.
        energies: Band energies in eV.
        degeneracies: Degeneracy of each band.
    """

    frac_coords: List[float]
    symbols: List[str]
    energies: List[float]
    degeneracies: List[int]

    @property
    def irrep_info_set(self) -> zip:
        """Iterate over (symbol, energy, degeneracy) tuples."""
        return zip(self.symbols, self.energies, self.degeneracies)


@dataclass
class Irreps(MSONable):
    """Collection of irreducible representations for multiple k-points.

    Attributes:
        sg_num: Space group number.
        irreps: Dictionary mapping k-point labels to Irrep objects.
                K-point for Gamma is "GM".
    """

    sg_num: int
    irreps: Dict[str, Irrep]

    def __call__(self) -> Dict[str, Irrep]:
        """Return the irreps dictionary."""
        return self.irreps

    def get_distances(self, x_ticks: "XTicks") -> List[List[float]]:
        """Get plot distances for each special point.

        Args:
            x_ticks: X-axis tick information.

        Returns:
            List of distance lists for each special point.
        """
        result: List[List[float]] = []
        for special_point in self.irreps:
            dist_list: List[float] = []
            for label, distance in zip(x_ticks.labels, x_ticks.distances):
                if special_point in label:
                    dist_list.append(distance)
            result.append(dist_list)
        return result


@dataclass(frozen=True)
class XTicks(MSONable):
    """X-axis tick information for band structure plots.

    Attributes:
        labels: Special point symbols (e.g., "Γ", "X", "C$\\mid$${\\rm C}_2$").
        distances: Corresponding cumulative distances along the path.
    """

    labels: List[str]
    distances: List[float]


@dataclass
class BandEdgeForPlot(MSONable):
    """Band edge information for plotting.

    Attributes:
        vbm: Valence band maximum energy.
        cbm: Conduction band minimum energy.
        vbm_distances: Distances where VBM occurs in the plot.
        cbm_distances: Distances where CBM occurs in the plot.
    """

    vbm: float
    cbm: float
    vbm_distances: List[float]
    cbm_distances: List[float]


# Backward compatibility alias
@dataclass
class BandEdge(BandEdgeForPlot):
    """Alias for BandEdgeForPlot (backward compatibility)."""
    pass


@dataclass
class BandEnergyInfo(MSONable):
    """Band energy data for a single calculation.

    Stores band energies organized by branch, spin, band index, and k-point.

    Attributes:
        band_energies: Nested list [branch][spin][band][kpoint].
                      Branches are separated by k-path discontinuities.
        band_edge: Band edge information (for semiconductors).
        fermi_level: Fermi level (for metals).
        irreps: Optional irreducible representations.
    """

    band_energies: List[List[List[List[float]]]]
    band_edge: Optional[BandEdgeForPlot] = None
    fermi_level: Optional[float] = None
    irreps: Optional[Irreps] = None

    def __post_init__(self) -> None:
        """Validate that either band_edge or fermi_level is provided."""
        if self.band_edge is None and self.fermi_level is None:
            raise ViseBandInfoError("Must provide either band_edge or fermi_level")

    def slide_energies(self, base_energy: float) -> None:
        """Shift all energies relative to base_energy."""
        self._slide_band_energies(base_energy)
        self._slide_band_edge(base_energy)
        self._slide_fermi_level(base_energy)
        self._slide_irreps(base_energy)

    def _slide_band_energies(self, base_energy: float) -> None:
        """Shift band energies."""
        self.band_energies = [
            [[[e - base_energy for e in kpts] for kpts in band]
             for band in spin]
            for spin in self.band_energies
        ]

    def _slide_band_edge(self, base_energy: float) -> None:
        """Shift band edge energies."""
        if self.band_edge:
            self.band_edge.vbm -= base_energy
            self.band_edge.cbm -= base_energy

    def _slide_fermi_level(self, base_energy: float) -> None:
        """Shift Fermi level."""
        if self.fermi_level:
            self.fermi_level -= base_energy

    def _slide_irreps(self, base_energy: float) -> None:
        """Shift irrep energies."""
        if self.irreps:
            for irrep in self.irreps.irreps.values():
                irrep.energies = [e - base_energy for e in irrep.energies]

    @property
    def is_magnetic(self) -> bool:
        """Whether calculation is spin-polarized."""
        return len(self.band_energies[0]) == 2

    def band_energy_region(
        self,
        decision_width: float = 0.1,
        bottom: Optional[float] = None,
        top: Optional[float] = None,
        offset: float = 0.0,
    ) -> List[List[float]]:
        """Estimate continuous band energy regions.

        Identifies energy ranges where bands form continuous groups,
        useful for determining valence/conduction band regions.

        Args:
            decision_width: Max energy gap to consider as continuous.
            bottom: Minimum energy to include.
            top: Maximum energy to include.
            offset: Value to subtract from all energies.

        Returns:
            List of [min, max] energy ranges.
        """
        result: List[List[float]] = []

        def add_boundary(lower: float, upper: float) -> None:
            result.append([lower - offset, upper - offset])

        # Flatten all band energies and sort
        sorted_energies = sorted(
            energy
            for branch in self.band_energies
            for spin in branch
            for band in spin
            for energy in band
        )

        # Apply energy filters
        if bottom is not None:
            sorted_energies = [e for e in sorted_energies if e >= bottom]
        if top is not None:
            sorted_energies = [e for e in sorted_energies if e <= top]

        prev_energy = sorted_energies.pop(0)
        lower_bound = prev_energy

        for energy in sorted_energies:
            if energy - prev_energy > decision_width:
                add_boundary(lower_bound, prev_energy)
                lower_bound = energy
            prev_energy = energy
        else:
            add_boundary(lower_bound, energy)

        return result


@dataclass
class BandPlotInfo(MSONable, ToJsonFileMixIn):
    """Complete band structure plot information.

    Supports multiple band structures with the same k-path for comparison.

    Attributes:
        band_energy_infos: Dictionary mapping names to BandEnergyInfo.
        distances_by_branch: Cumulative distances for each k-path branch.
        x_ticks: X-axis tick information.
        title: Overall plot title.
    """

    band_energy_infos: Dict[str, BandEnergyInfo]
    distances_by_branch: List[List[float]]
    x_ticks: XTicks
    title: Optional[str] = None

    def __post_init__(self) -> None:
        """Validate distance consistency."""
        assert self.distances_by_branch[0][0] == self.x_ticks.distances[0]
        assert self.distances_by_branch[-1][-1] == self.x_ticks.distances[-1]

    def __add__(self, other: "BandPlotInfo") -> "BandPlotInfo":
        """Merge two BandPlotInfo objects."""
        assert self.distances_by_branch == other.distances_by_branch
        merged_infos = deepcopy(self.band_energy_infos)
        merged_infos.update(other.band_energy_infos)
        return BandPlotInfo(
            merged_infos, self.distances_by_branch, self.x_ticks, self.title
        )


class ViseBandInfoError(ViseError):
    """Error raised for invalid band information."""
    pass


class BandMplSettings:
    """Matplotlib styling settings for band structure plots.

    Attributes:
        colors: Color cycle for different band structures.
        linewidth: Line width cycle.
        circle_size: Size of band edge markers.
        hline: Settings for horizontal lines (band edges).
        title_font_size: Font size for title.
        label_font_size: Font size for axis labels.
        show_legend: Whether to display legend.
        legend: Legend configuration.
    """

    def __init__(
        self,
        colors: Optional[List[str]] = None,
        linewidth: Optional[List[float]] = None,
        circle_size: int = 70,
        circle_colors: Optional[List[str]] = None,
        band_edge_line_width: float = 0.75,
        band_edge_line_color: str = "black",
        band_edge_line_style: str = "-.",
        title_font_size: int = 15,
        label_font_size: int = 15,
        show_legend: bool = True,
        legend_location: str = "lower right",
    ) -> None:
        """Initialize band plot settings."""
        self.colors = colors or ["#E15759", "#4E79A7", "#F28E2B", "#76B7B2"]
        self.linewidth = cycle(linewidth) if linewidth else cycle([1.0])

        self.circle_size = circle_size
        self.hline = {
            "linewidth": band_edge_line_width,
            "color": band_edge_line_color,
            "linestyle": band_edge_line_style,
        }

        self.title_font_size = title_font_size
        self.label_font_size = label_font_size

        self.circle_colors = circle_colors or ["pink", "green"]

        self.show_legend = show_legend
        self.legend = {"loc": legend_location}

    def band_structure(self, index: int, label: str) -> dict:
        """Get line style settings for a band structure."""
        return {
            "color": self.colors[index],
            "linewidth": next(self.linewidth),
            "label": label,
        }

    def circle(self, index: int) -> dict:
        """Get marker settings for band edge points."""
        return {
            "color": self.colors[index],
            "marker": "o",
            "s": self.circle_size,
        }


class BandMplPlotter:
    """Matplotlib-based band structure plotter.

    Creates band structure plots with support for multiple calculations,
    band edge markers, and irreducible representations.

    Attributes:
        plt: Reference to matplotlib.pyplot.
        title: Plot title.
        mpl_defaults: Style settings.
    """

    def __init__(
        self,
        band_plot_info: BandPlotInfo,
        energy_range: List[float],
        base_energy: Optional[float] = None,
        base_energy_title: Optional[str] = None,
        mpl_defaults: Optional[BandMplSettings] = None,
    ) -> None:
        """Initialize the band plotter.

        Args:
            band_plot_info: Band structure data.
            energy_range: Y-axis energy range [min, max].
            base_energy: Reference energy for alignment.
            base_energy_title: Which calculation to use for base energy.
            mpl_defaults: Matplotlib styling settings.
        """
        # Deep copy to avoid side effects from energy shifts
        band_plot_info = deepcopy(band_plot_info)
        self.band_energy_infos = band_plot_info.band_energy_infos
        self.distances_by_branch = band_plot_info.distances_by_branch
        self.x_ticks = band_plot_info.x_ticks

        self.energy_range = energy_range
        self.title = band_plot_info.title
        self.mpl_defaults = mpl_defaults or BandMplSettings()
        self.plt = plt

        if base_energy is None:
            base_energy = get_base_energy(self.band_energy_infos, base_energy_title)
        slide_band_energies(self.band_energy_infos, base_energy)

    def construct_plot(self) -> None:
        """Build the complete band structure plot."""
        self._add_band_set()
        self._set_figure_legend()
        self._set_x_range()
        self._set_y_range()
        self._set_labels()
        self._set_x_ticks()
        self._set_title()
        self._set_formatter()
        self.plt.tight_layout()

    def _add_band_set(self) -> None:
        """Add all band structures to the plot."""
        for index, (name, info) in enumerate(self.band_energy_infos.items()):
            self._add_band_structures(info, name, index)

            if info.band_edge:
                self._add_band_edge(info.band_edge, index)
            elif info.fermi_level:
                self._add_fermi_level(info.fermi_level)

    def _add_band_structures(
        self, band_info: BandEnergyInfo, band_name: str, index: int
    ) -> None:
        """Add band structure curves."""
        mpl_args = self.mpl_defaults.band_structure(index, band_name)

        for distances, branch_energies in zip(
            self.distances_by_branch, band_info.band_energies
        ):
            for spin_idx, spin_energies in enumerate(branch_energies):
                # Use dotted line for minority spin
                if spin_idx == 0:
                    mpl_args.pop("linestyle", None)
                else:
                    mpl_args["linestyle"] = ":"

                for band_energies in spin_energies:
                    self.plt.plot(distances, band_energies, **mpl_args)
                    mpl_args.pop("label", None)  # Only label first curve

    def _add_band_edge(self, band_edge: BandEdgeForPlot, index: int) -> None:
        """Add VBM/CBM lines and markers."""
        self.plt.axhline(y=band_edge.vbm, **self.mpl_defaults.hline)
        self.plt.axhline(y=band_edge.cbm, **self.mpl_defaults.hline)

        for dist in band_edge.vbm_distances:
            self.plt.scatter(
                dist, band_edge.vbm, **self.mpl_defaults.circle(index)
            )
        for dist in band_edge.cbm_distances:
            self.plt.scatter(
                dist, band_edge.cbm, **self.mpl_defaults.circle(index)
            )

    def _add_fermi_level(self, fermi_level: float) -> None:
        """Add Fermi level line for metals."""
        self.plt.axhline(y=fermi_level, **self.mpl_defaults.hline)

    def _set_figure_legend(self) -> None:
        """Display legend if multiple band structures."""
        if self.mpl_defaults.show_legend and len(self.band_energy_infos) > 1:
            self.plt.legend(**self.mpl_defaults.legend)

    def _set_x_range(self) -> None:
        """Set x-axis range."""
        self.plt.xlim(self.x_ticks.distances[0], self.x_ticks.distances[-1])

    def _set_y_range(self) -> None:
        """Set y-axis energy range."""
        self.plt.ylim(self.energy_range[0], self.energy_range[1])

    def _set_labels(self) -> None:
        """Set axis labels."""
        self.plt.xlabel("Wave vector", size=self.mpl_defaults.label_font_size)
        self.plt.ylabel("Energy (eV)", size=self.mpl_defaults.label_font_size)

    def _set_x_ticks(self) -> None:
        """Set x-axis ticks and vertical lines."""
        axis = self.plt.gca()
        axis.set_xticks(self.x_ticks.distances)
        axis.set_xticklabels(self.x_ticks.labels)

        for distance, label in zip(
            self.x_ticks.distances[1:-1], self.x_ticks.labels[1:-1]
        ):
            linestyle = "-" if "\\mid" in label else "--"
            plt.axvline(x=distance, linestyle=linestyle)

    def _set_title(self) -> None:
        """Set plot title."""
        if self.title:
            self.plt.title(self.title, size=self.mpl_defaults.title_font_size)

    def _set_formatter(self) -> None:
        """Set axis formatters."""
        axis = self.plt.gca()
        axis.yaxis.set_major_formatter(float_to_int_formatter)


def get_base_energy(
    band_energy_infos: Dict[str, BandEnergyInfo],
    base_energy_title: Optional[str] = None,
) -> float:
    """Determine base energy for alignment.

    Args:
        band_energy_infos: Dictionary of band structure data.
        base_energy_title: Which calculation to use (None for first).

    Returns:
        Base energy (VBM for semiconductors, Fermi level for metals).
    """
    if base_energy_title is None:
        base_energy_title = next(iter(band_energy_infos))
        if len(band_energy_infos) > 1:
            logger.warning(f"Base energy set from '{base_energy_title}'.")

    base_info = band_energy_infos[base_energy_title]

    if base_info.band_edge:
        return base_info.band_edge.vbm
    return base_info.fermi_level


def slide_band_energies(
    band_energy_infos: Dict[str, BandEnergyInfo],
    base_energy: float,
) -> None:
    """Shift all band energies relative to base.

    Args:
        band_energy_infos: Dictionary of band structure data.
        base_energy: Energy to use as zero reference.
    """
    for band_info in band_energy_infos.values():
        band_info.slide_energies(base_energy)
