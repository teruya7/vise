# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Dielectric function analysis API.

This module provides functions for analyzing dielectric functions
from VASP output files.

Example:
    >>> from vise.api.dielectric import analyze_dielectric
    >>>
    >>> result = analyze_dielectric("vasprun.xml", "OUTCAR")
    >>> result.save_plot("diele.pdf")
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple, Union

from vise.analyzer.dielectric_function import DieleFuncData
from vise.analyzer.vasp.make_diele_func import make_diele_func
from vise.analyzer.plot_diele_func_data import DieleFuncMplPlotter, DieleFuncPlotType


@dataclass
class DielectricAnalysisResult:
    """Dielectric function analysis result.
    
    Attributes:
        diele_func_data: DieleFuncData object containing dielectric data
    """
    diele_func_data: DieleFuncData
    
    def save_plot(
        self,
        filename: Union[str, Path] = "dielectric.pdf",
        plot_type: Union[str, DieleFuncPlotType] = "real",
        directions: Optional[List[str]] = None,
        y_range: Optional[Tuple[float, float]] = None,
        title: Optional[str] = None,
        format: str = "pdf"
    ) -> None:
        """
        Save dielectric function plot to file.
        
        Args:
            filename: Output filename
            plot_type: Type of plot:
                - "real": Real part of dielectric function
                - "imag": Imaginary part
                - "absorption_coeff": Absorption coefficient
            directions: Directions to plot (e.g., ["xx", "yy", "zz"])
            y_range: Y-axis range
            title: Plot title
            format: Output format (pdf, png, etc.)
        """
        if isinstance(plot_type, str):
            plot_type = DieleFuncPlotType.from_string(plot_type)
        
        y_range_processed = None
        if y_range:
            if plot_type is DieleFuncPlotType.absorption_coeff:
                y_range_processed = [10 ** y_range[0], 10 ** y_range[1]]
            else:
                y_range_processed = list(y_range)
        
        if title:
            self.diele_func_data.title = title
        
        plotter = DieleFuncMplPlotter(self.diele_func_data)
        plotter.construct_plot(
            directions=directions,
            plot_type=plot_type,
            y_range=y_range_processed
        )
        plotter.plt.savefig(str(filename), format=format)
        plotter.plt.close()
    
    def save_json(self, filename: Union[str, Path] = "diele_func_data.json") -> None:
        """Save dielectric function data as JSON file."""
        self.diele_func_data.to_json_file(str(filename))
    
    def save_csv(self, filename: Union[str, Path] = "diele_func_data.csv") -> None:
        """Save dielectric function data as CSV file."""
        self.diele_func_data.to_csv_file(str(filename))


def analyze_dielectric(
    vasprun: Union[str, Path],
    outcar: Union[str, Path] = "OUTCAR",
    use_vasp_real: bool = True,
    ita: float = 0.01
) -> DielectricAnalysisResult:
    """
    Analyze dielectric function from VASP output.
    
    Args:
        vasprun: Path to vasprun.xml file
        outcar: Path to OUTCAR file
        use_vasp_real: Use real part from VASP (True) or calculate via
            Kramers-Kronig (False)
        ita: Broadening parameter for Kramers-Kronig transformation
    
    Returns:
        DielectricAnalysisResult with dielectric function data
    
    Example:
        >>> from vise.api.dielectric import analyze_dielectric
        >>>
        >>> result = analyze_dielectric("vasprun.xml", "OUTCAR")
        >>> result.save_plot("real.pdf", plot_type="real")
        >>> result.save_plot("imag.pdf", plot_type="imag")
        >>> result.save_plot("absorption.pdf", plot_type="absorption_coeff")
    """
    from pymatgen.io.vasp import Vasprun, Outcar
    
    vasprun_obj = Vasprun(str(vasprun))
    outcar_obj = Outcar(str(outcar))
    
    diele_func_data = make_diele_func(
        vasprun_obj,
        outcar_obj,
        use_vasp_real=use_vasp_real,
        ita=ita
    )
    
    return DielectricAnalysisResult(diele_func_data=diele_func_data)


def load_dielectric_from_csv(csv_filename: Union[str, Path]) -> DielectricAnalysisResult:
    """
    Load dielectric function data from CSV file.
    
    Args:
        csv_filename: Path to CSV file
    
    Returns:
        DielectricAnalysisResult with loaded data
    """
    diele_func_data = DieleFuncData.from_csv_file(str(csv_filename))
    return DielectricAnalysisResult(diele_func_data=diele_func_data)
