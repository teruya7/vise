# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Band structure analysis API.

This module provides functions for analyzing band structure
from VASP output files.

Example:
    >>> from vise.api.band import analyze_band
    >>>
    >>> result = analyze_band("vasprun.xml", "KPOINTS")
    >>> print(f"Band gap: {result.band_gap} eV")
    >>> result.save_plot("band.pdf")
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Union

from vise.analyzer.plot_band import BandMplPlotter, BandPlotInfo
from vise.analyzer.vasp.plot_band import BandPlotInfoFromVasp


@dataclass
class BandAnalysisResult:
    """Band structure analysis result.
    
    Attributes:
        plot_info: BandPlotInfo object containing band data
        band_gap: Band gap in eV (None if metallic)
        vbm: Valence band maximum energy in eV
        cbm: Conduction band minimum energy in eV
        is_metal: True if no band gap exists
    """
    plot_info: BandPlotInfo
    band_gap: Optional[float]
    vbm: Optional[float]
    cbm: Optional[float]
    is_metal: bool

    def save_plot(
        self,
        filename: Union[str, Path] = "band.pdf",
        energy_range: Tuple[float, float] = (-10.0, 10.0),
        format: str = "pdf"
    ) -> None:
        """
        Save band structure plot to file.
        
        Args:
            filename: Output filename
            energy_range: Energy range (min, max) for y-axis
            format: Output format (pdf, png, etc.)
        """
        plotter = BandMplPlotter(self.plot_info, energy_range=list(energy_range))
        plotter.construct_plot()
        plotter.plt.savefig(str(filename), format=format)
        plotter.plt.close()

    def save_json(self, filename: Union[str, Path] = "band_plot_info.json") -> None:
        """
        Save plot info as JSON file.
        
        Args:
            filename: Output JSON filename
        """
        self.plot_info.to_json_file(str(filename))


def analyze_band(
    vasprun: Union[str, Path],
    kpoints_filename: Union[str, Path] = "KPOINTS"
) -> BandAnalysisResult:
    """
    Analyze band structure from VASP output.
    
    Args:
        vasprun: Path to vasprun.xml file
        kpoints_filename: Path to KPOINTS file
    
    Returns:
        BandAnalysisResult with band structure information
    
    Example:
        >>> from vise.api.band import analyze_band
        >>>
        >>> result = analyze_band("vasprun.xml")
        >>> if result.is_metal:
        ...     print("System is metallic")
        ... else:
        ...     print(f"Band gap: {result.band_gap:.2f} eV")
        >>> result.save_plot("my_band.pdf", energy_range=(-5, 5))
    """
    from pymatgen.io.vasp import Vasprun

    vasprun_obj = Vasprun(str(vasprun))

    band_plot_info_from_vasp = BandPlotInfoFromVasp(
        vasprun=vasprun_obj,
        kpoints_filename=str(kpoints_filename)
    )
    plot_info = band_plot_info_from_vasp.make_band_plot_info()

    band_edge = plot_info.band_edge

    if band_edge is None:
        return BandAnalysisResult(
            plot_info=plot_info,
            band_gap=None,
            vbm=None,
            cbm=None,
            is_metal=True
        )

    return BandAnalysisResult(
        plot_info=plot_info,
        band_gap=band_edge.band_gap if hasattr(band_edge, 'band_gap') else None,
        vbm=band_edge.vbm_info.energy if band_edge.vbm_info else None,
        cbm=band_edge.cbm_info.energy if band_edge.cbm_info else None,
        is_metal=band_edge.band_gap is None if hasattr(band_edge, 'band_gap') else True
    )
