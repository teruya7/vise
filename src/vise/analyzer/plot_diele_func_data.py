# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Dielectric function plotting utilities.

This module provides classes for plotting dielectric function data
including absorption coefficient, refractive index, and reflectivity,
using matplotlib as the backend.
"""

from abc import abstractmethod
from typing import List, Optional, Tuple

import numpy as np
from matplotlib import pyplot as plt

from vise.analyzer.dielectric_function import DieleFuncData
from vise.util.enum import ExtendedEnum
from vise.util.matplotlib import float_to_int_formatter


class TensorDirection(ExtendedEnum):
    """Tensor direction for extracting specific components.

    For 3x3 dielectric tensors, select average or specific components.
    """

    average = "average"
    xx = "xx"
    yy = "yy"
    zz = "zz"
    xy = "xy"
    yz = "yz"
    xz = "xz"

    def val(self, tensors: List[List[float]]) -> List[float]:
        """Extract values for this direction from tensor data.

        Args:
            tensors: List of 6-component tensors [xx, yy, zz, xy, yz, xz].

        Returns:
            List of values for the specified direction.
        """
        if self is self.average:
            return [(t[0] + t[1] + t[2]) / 3 for t in tensors]
        else:
            component_order = ["xx", "yy", "zz", "xy", "yz", "xz"]
            idx = component_order.index(self.name)
            return [tensor[idx] for tensor in tensors]


class DieleFuncPlotType(ExtendedEnum):
    """Type of dielectric function plot.

    Determines what physical quantity to plot and how to format axes.
    """

    diele_func = "diele_func"
    absorption_coeff = "absorption_coeff"
    refraction = "refraction"
    reflectivity = "reflectivity"

    def y_axis_label(self, plot_engine: str) -> str:
        """Get formatted y-axis label for this plot type.

        Args:
            plot_engine: 'plotly' or 'matplotlib' for appropriate formatting.

        Returns:
            Y-axis label with units if applicable.

        Raises:
            ValueError: If plot type or engine not recognized.
        """
        labels = {
            self.diele_func: "Dielectric function",
            self.absorption_coeff: "Absorption coefficient",
            self.refraction: "Refraction",
            self.reflectivity: "Reflectivity",
        }

        label = labels.get(self)
        if label is None:
            raise ValueError(f"Label not implemented for {self}")

        # Add units for absorption coefficient
        if self is self.absorption_coeff:
            if plot_engine == "plotly":
                return f"{label} (cm<sup>-1</sup>)"
            elif plot_engine == "matplotlib":
                return f"{label} (cm$^{{-1}}$)"
            else:
                raise ValueError(f"Unknown plot engine: {plot_engine}")

        return label

    @property
    def y_axis_default_limits(self) -> Optional[Tuple[float, float]]:
        """Default y-axis limits for this plot type."""
        if self is self.absorption_coeff:
            return (10**3, 10**7)
        return None

    def tensors(self, data: DieleFuncData) -> Tuple[List, ...]:
        """Extract relevant tensor data from DieleFuncData.

        Args:
            data: Dielectric function data object.

        Returns:
            Tuple of tensor lists (real, imag) or single tensor list.
        """
        if self is self.diele_func:
            return data.diele_func_real, data.diele_func_imag
        elif self is self.absorption_coeff:
            return (data.absorption_coeff,)
        elif self is self.refraction:
            return data.refractive_idx_real, data.refractive_idx_imag
        elif self is self.reflectivity:
            return (data.reflectivity,)


def auto_y_range(
    plot_type: DieleFuncPlotType,
    values: List[float],
) -> Tuple[float, float]:
    """Calculate automatic y-axis range for plot values.

    Args:
        plot_type: Type of dielectric function plot.
        values: All values to be plotted.

    Returns:
        Tuple of (y_min, y_max).
    """
    max_val = max(values)

    if plot_type is DieleFuncPlotType.absorption_coeff:
        return (10**3, max_val * 2)

    min_val = min(values)
    # Add 5% padding on each side
    return (1.05 * min_val - 0.05 * max_val, 1.05 * max_val - 0.05 * min_val)


class DieleFuncPlotter:
    """Base class for dielectric function plotting.

    Provides common functionality for different plotting backends.
    """

    def __init__(
        self,
        diele_func_data: DieleFuncData,
        energy_range: Optional[List[float]] = None,
    ) -> None:
        """Initialize the plotter.

        Args:
            diele_func_data: Dielectric function data to plot.
            energy_range: Energy range [min, max] in eV.
        """
        self.diele_func_data = diele_func_data
        self.energies = diele_func_data.energies
        self.band_gap = diele_func_data.band_gap
        self.energy_range = energy_range or [0, 10]
        self._x_axis_title = "Energy (eV)"

    def add_plot(
        self,
        directions: List[str],
        plot_type: DieleFuncPlotType,
    ) -> List[float]:
        """Add plots for specified directions and collect all values.

        Args:
            directions: List of direction labels to plot.
            plot_type: Type of plot to generate.

        Returns:
            All plotted values (for y-range calculation).
        """
        max_energy_idx = np.where(
            np.array(self.energies) > self.energy_range[1]
        )[0][0]

        all_values: List[float] = []
        tensor_data = plot_type.tensors(self.diele_func_data)

        for component_idx, (component_name, tensors) in enumerate(
            zip(["real", "imag"], tensor_data)
        ):
            for direction in directions:
                dir_idx = self.diele_func_data.directions.index(direction)
                tensor = tensors[dir_idx]
                all_values.extend(tensor[:max_energy_idx])

                # Name: direction only if single component, else with real/imag
                if len(tensor_data) == 1:
                    name = direction
                else:
                    name = f"{component_name}_{direction}"

                self.add_single_plot(name, tensor, plot_type)

        return all_values

    @abstractmethod
    def add_single_plot(
        self,
        name: str,
        tensor: List[float],
        plot_type: DieleFuncPlotType,
    ) -> None:
        """Add a single curve to the plot. Override in subclasses."""
        pass


class DieleFuncMplPlotter(DieleFuncPlotter):
    """Matplotlib-based dielectric function plotter."""

    def __init__(
        self,
        diele_func_data: DieleFuncData,
        energy_range: Optional[List[float]] = None,
    ) -> None:
        """Initialize matplotlib plotter."""
        super().__init__(diele_func_data, energy_range)
        self.title = diele_func_data.title
        self.plt = plt
        self.plt.clf()

    def construct_plot(
        self,
        directions: Tuple[str, ...] = ("ave",),
        plot_type: DieleFuncPlotType = DieleFuncPlotType.absorption_coeff,
        y_range: Optional[Tuple[float, float]] = None,
    ) -> None:
        """Build the complete plot.

        Args:
            directions: Tensor directions to plot.
            plot_type: Type of dielectric function plot.
            y_range: Optional explicit y-axis range.
        """
        self._add_curves(directions, plot_type, y_range)
        self._add_band_gap_line()
        self._set_legend()
        self._set_x_range()
        self._set_axis_labels(plot_type)
        self._set_formatter()
        self._set_title()
        self.plt.tight_layout()

    def _add_curves(
        self,
        directions: Tuple[str, ...],
        plot_type: DieleFuncPlotType,
        y_range: Optional[Tuple[float, float]],
    ) -> None:
        """Add all data curves and set y-range."""
        all_values = self.add_plot(list(directions), plot_type)
        y_min, y_max = y_range or auto_y_range(plot_type, all_values)
        self.plt.gca().set_ylim(ymin=y_min, ymax=y_max)

    def add_single_plot(
        self,
        name: str,
        tensor: List[float],
        plot_type: DieleFuncPlotType,
    ) -> None:
        """Add a single curve to the matplotlib plot."""
        if plot_type is DieleFuncPlotType.absorption_coeff:
            self.plt.semilogy(self.energies, tensor, label=name)
        else:
            self.plt.plot(self.energies, tensor, label=name)

    def _add_band_gap_line(self) -> None:
        """Add vertical line at band gap energy."""
        self.plt.axvline(
            x=self.band_gap, linestyle="dashed", color="black", linewidth=1
        )

    def _set_legend(self) -> None:
        """Display figure legend."""
        self.plt.legend()

    def _set_x_range(self) -> None:
        """Set x-axis energy range."""
        self.plt.xlim(self.energy_range[0], self.energy_range[1])

    def _set_axis_labels(self, plot_type: DieleFuncPlotType) -> None:
        """Set axis labels."""
        self.plt.xlabel(self._x_axis_title)
        self.plt.ylabel(plot_type.y_axis_label("matplotlib"))

    def _set_formatter(self) -> None:
        """Set axis tick formatters."""
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)

    def _set_title(self) -> None:
        """Set plot title if available."""
        if self.title:
            self.plt.title(self.title)
