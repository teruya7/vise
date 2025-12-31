# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Density of States (DOS) analysis API.

This module provides functions for analyzing DOS from VASP output files.

Example:
    >>> from vise.api.dos import analyze_dos
    >>>
    >>> result = analyze_dos("vasprun.xml", "OUTCAR")
    >>> print(f"Band gap: {result.band_gap} eV")
    >>> result.save_plot("dos.pdf")
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple, Union, Dict

from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.analyzer.vasp.dos_data import DosDataFromVasp
from vise.analyzer.dos_data import DosData
from vise.analyzer.plot_dos import DosPlotter, DosPlotData
from vise.analyzer.atom_grouping_type import AtomGroupingType


@dataclass
class DosAnalysisResult:
    """DOS analysis result.
    
    Attributes:
        dos_data: DosData object containing DOS information
        band_gap: Band gap in eV (None if metallic)
        vbm: Valence band maximum energy in eV
        cbm: Conduction band minimum energy in eV  
        is_metal: True if no band gap exists
        efermi: Fermi energy in eV
    """
    dos_data: DosData
    band_gap: Optional[float]
    vbm: Optional[float]
    cbm: Optional[float]
    is_metal: bool
    efermi: float
    _structure: any = None
    _vertical_lines: List[float] = None
    
    def save_plot(
        self,
        filename: Union[str, Path] = "dos.pdf",
        energy_range: Tuple[float, float] = (-10.0, 10.0),
        grouping_type: Union[str, AtomGroupingType] = "atoms",
        target: Optional[List[str]] = None,
        y_max_ranges: Optional[List[float]] = None,
        title: Optional[str] = None,
        show_legend: bool = True,
        format: str = "pdf"
    ) -> None:
        """
        Save DOS plot to file.
        
        Args:
            filename: Output filename
            energy_range: Energy range (min, max) for x-axis
            grouping_type: How to group atoms ("atoms", "elements", etc.)
            target: Target atoms/elements to plot (None for all)
            y_max_ranges: Maximum y-axis values per subplot
            title: Plot title
            show_legend: Whether to show legend
            format: Output format (pdf, png, etc.)
        """
        if isinstance(grouping_type, str):
            grouping_type = AtomGroupingType.from_string(grouping_type)
        
        grouped_atom_indices = grouping_type.grouped_atom_indices(
            self._structure, target
        )
        
        ylim_set = None
        if y_max_ranges:
            if self.dos_data.spin:
                ylim_set = [[-y_max, y_max] for y_max in y_max_ranges]
            else:
                ylim_set = [[0, y_max] for y_max in y_max_ranges]
        
        plot_data = self.dos_data.dos_plot_data(
            grouped_atom_indices,
            energy_range=list(energy_range),
            dos_ranges=ylim_set,
            title=title
        )
        
        plotter = DosPlotter(plot_data, show_legend)
        plotter.construct_plot()
        plotter.plt.savefig(str(filename), format=format)
        plotter.plt.close()
    
    def save_json(self, filename: Union[str, Path] = "dos_data.json") -> None:
        """Save DOS data as JSON file."""
        self.dos_data.to_json_file(str(filename))


def analyze_dos(
    vasprun: Union[str, Path],
    outcar: Union[str, Path] = "OUTCAR",
    base_energy: Optional[float] = None,
    crop_first_value: bool = True
) -> DosAnalysisResult:
    """
    Analyze DOS from VASP output.
    
    Args:
        vasprun: Path to vasprun.xml file
        outcar: Path to OUTCAR file
        base_energy: Reference energy for alignment (default: VBM or Fermi)
        crop_first_value: Whether to crop the first DOS value
    
    Returns:
        DosAnalysisResult with DOS information
    
    Example:
        >>> from vise.api.dos import analyze_dos
        >>>
        >>> result = analyze_dos("vasprun.xml", "OUTCAR")
        >>> print(f"Fermi energy: {result.efermi:.2f} eV")
        >>> result.save_plot("dos.pdf", energy_range=(-5, 5))
    """
    from pymatgen.io.vasp import Vasprun, Outcar
    
    vasprun_obj = Vasprun(str(vasprun))
    outcar_obj = Outcar(str(outcar))
    
    band_edge = VaspBandEdgeProperties(vasprun_obj, outcar_obj)
    
    if band_edge.band_gap:
        vertical_lines = [band_edge.vbm_info.energy, band_edge.cbm_info.energy]
        vbm = band_edge.vbm_info.energy
        cbm = band_edge.cbm_info.energy
        band_gap = band_edge.band_gap
        is_metal = False
    else:
        vertical_lines = [vasprun_obj.efermi]
        vbm = None
        cbm = None
        band_gap = None
        is_metal = True
    
    if base_energy is None:
        base = vertical_lines[0]
    else:
        base = base_energy
    
    dos_data_from_vasp = DosDataFromVasp(
        vasprun_obj, vertical_lines, base, crop_first_value
    )
    dos_data = dos_data_from_vasp.make_dos_data()
    
    result = DosAnalysisResult(
        dos_data=dos_data,
        band_gap=band_gap,
        vbm=vbm,
        cbm=cbm,
        is_metal=is_metal,
        efermi=vasprun_obj.efermi
    )
    result._structure = vasprun_obj.final_structure
    result._vertical_lines = vertical_lines
    
    return result
