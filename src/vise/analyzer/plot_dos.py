# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Density of states (DOS) plotting utilities.

This module provides matplotlib-based plotting for electronic density
of states, supporting multiple subplot panels, spin-polarized data,
and customizable styling.
"""

from typing import List, Optional

import matplotlib.pyplot as plt

from vise.analyzer.dos_data import DosPlotData
from vise.util.matplotlib import float_to_int_formatter


class DosMplSettings:
    """Matplotlib styling settings for DOS plots.

    Attributes:
        colors: Color cycle for DOS lines.
        linewidth: Line width for DOS curves.
        vline: Settings dict for vertical band edge lines.
        title_font_size: Font size for plot titles.
        label_font_size: Font size for axis labels.
    """

    def __init__(
        self,
        colors: Optional[List[str]] = None,
        linewidth: float = 1.0,
        band_edge_line_width: float = 0.75,
        band_edge_line_color: str = "black",
        band_edge_line_style: str = "-.",
        title_font_size: int = 15,
        label_font_size: int = 12,
    ) -> None:
        """Initialize DOS plot settings.

        Args:
            colors: List of colors for DOS lines. Defaults to a 4-color palette.
            linewidth: Width of DOS lines.
            band_edge_line_width: Width of VBM/CBM marker lines.
            band_edge_line_color: Color of band edge lines.
            band_edge_line_style: Line style for band edges.
            title_font_size: Font size for title.
            label_font_size: Font size for axis labels.
        """
        self.colors = colors or ["#36454f", "#E15759", "#4E79A7", "#F28E2B"]
        self.linewidth = linewidth
        self.vline = {
            "linewidth": band_edge_line_width,
            "color": band_edge_line_color,
            "linestyle": band_edge_line_style,
        }
        self.title_font_size = title_font_size
        self.label_font_size = label_font_size

    def dos_line(self, index: int) -> dict:
        """Get line style settings for DOS at given index."""
        return {
            "color": self.colors[index % len(self.colors)],
            "linewidth": self.linewidth,
        }


class DosPlotter:
    """Matplotlib-based DOS plotter.

    Creates multi-panel DOS plots with support for spin-polarized data,
    legends, and band edge markers.

    Attributes:
        fig: Matplotlib figure object.
        plt: Reference to matplotlib.pyplot.
        title: Plot title from data.
        mpl_defaults: Style settings.
    """

    def __init__(
        self,
        dos_plot_data: DosPlotData,
        show_legend: bool = True,
        mpl_defaults: Optional[DosMplSettings] = None,
    ) -> None:
        """Initialize the DOS plotter.

        Args:
            dos_plot_data: DOS data to plot.
            show_legend: Whether to show legends on plots.
            mpl_defaults: Matplotlib styling settings.
        """
        self._dos_plot_data = dos_plot_data
        self._show_legend = show_legend
        self.mpl_defaults = mpl_defaults or DosMplSettings()
        self.title = dos_plot_data.title
        self.plt = plt

        num_panels = len(self._dos_plot_data.doses)
        self.fig, self._axs = self.plt.subplots(
            num_panels, 1, sharex=True, gridspec_kw={"hspace": 0.1}
        )

        # Ensure _axs is always a list
        if num_panels == 1:
            self._axs = [self._axs]

    def construct_plot(self) -> None:
        """Build the complete DOS plot."""
        self._axs[0].set_xlim(self._dos_plot_data.energy_range)
        self._axs[0].xaxis.set_major_formatter(float_to_int_formatter)

        for i in range(len(self._dos_plot_data.doses)):
            self._add_panel(i)

        self._set_x_labels()
        middle_idx = int(i / 2)
        self._set_ylabel(middle_idx)
        self._set_title()

    def _add_panel(self, panel_idx: int) -> None:
        """Configure a single DOS panel."""
        self._add_dos_curves(panel_idx)
        self._set_y_range(panel_idx)
        self._set_dos_zero_line(panel_idx)
        self._set_band_edge_lines(panel_idx)
        self._set_formatter(panel_idx)

    def _add_dos_curves(self, panel_idx: int) -> None:
        """Add DOS curves to a panel."""
        panel_data = self._dos_plot_data.doses[panel_idx]

        for dos_idx, named_dos in enumerate(panel_data):
            for spin_idx, spin_dos in enumerate(named_dos.dos):
                # Up spin: positive, Down spin: negative
                sign = 1 - 2 * spin_idx
                dos_values = [d * sign for d in spin_dos]

                line_kwargs = self.mpl_defaults.dos_line(dos_idx)

                # Add label only for first spin component
                if spin_idx == 0 and self._show_legend:
                    if named_dos.name:
                        label = f"{self._dos_plot_data.names[panel_idx]}-{named_dos.name}"
                    else:
                        label = self._dos_plot_data.names[panel_idx]
                    line_kwargs["label"] = label

                self._axs[panel_idx].plot(
                    self._dos_plot_data.relative_energies,
                    dos_values,
                    **line_kwargs,
                )

        if self._show_legend:
            self._axs[panel_idx].legend(
                bbox_to_anchor=(0.9, 1),
                loc="upper left",
                borderaxespad=0,
                markerscale=0.1,
            )

    def _set_ylabel(self, panel_idx: int) -> None:
        """Set y-axis label on the specified panel."""
        self._axs[panel_idx].set_ylabel(
            "Dos (1/eV)", size=self.mpl_defaults.label_font_size
        )

    def _set_dos_zero_line(self, panel_idx: int) -> None:
        """Add horizontal line at DOS=0."""
        self._axs[panel_idx].axhline(0, linestyle=":", color="black")

    def _set_y_range(self, panel_idx: int) -> None:
        """Set y-axis range for a panel."""
        self._axs[panel_idx].set_ylim(self._dos_plot_data.dos_ranges[panel_idx])

    def _set_band_edge_lines(self, panel_idx: int) -> None:
        """Add vertical lines at VBM/CBM positions."""
        for energy in self._dos_plot_data.energy_lines:
            self._axs[panel_idx].axvline(x=energy, **self.mpl_defaults.vline)

    def _set_formatter(self, panel_idx: int) -> None:
        """Set axis formatters."""
        self._axs[panel_idx].yaxis.set_major_formatter(float_to_int_formatter)

    def _set_x_labels(self) -> None:
        """Set x-axis label on the bottom panel."""
        self.plt.xlabel("Energy (eV)", size=self.mpl_defaults.label_font_size)

    def _set_title(self) -> None:
        """Set figure title if available."""
        if self.title:
            if self.fig:
                self.fig.suptitle(
                    self.title, fontsize=self.mpl_defaults.title_font_size
                )
                self.fig.subplots_adjust(top=0.92)
            else:
                self.plt.title(
                    self.title, fontsize=self.mpl_defaults.title_font_size
                )
