# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Band edge and effective mass analysis API.

This module provides functions for analyzing band edge properties
and effective mass from VASP output files.

Example:
    >>> from vise.api.band_edge import get_band_edge_properties
    >>>
    >>> props = get_band_edge_properties("vasprun.xml", "OUTCAR")
    >>> print(f"Band gap: {props.band_gap} eV")
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple, Union

from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.analyzer.vasp.make_effective_mass import make_effective_mass


@dataclass
class KpointInfo:
    """K-point information for band edge.
    
    Attributes:
        kpoint_index: Index of k-point
        kpoint_coords: Fractional coordinates of k-point
        energy: Energy at this k-point in eV
    """
    kpoint_index: int
    kpoint_coords: Tuple[float, float, float]
    energy: float


@dataclass
class BandEdgeResult:
    """Band edge properties result.
    
    Attributes:
        band_gap: Band gap in eV (None if metallic)
        vbm_info: VBM k-point information
        cbm_info: CBM k-point information
        is_direct: True if direct band gap
        is_metal: True if no band gap
        efermi: Fermi energy in eV
    """
    band_gap: Optional[float]
    vbm_info: Optional[KpointInfo]
    cbm_info: Optional[KpointInfo]
    is_direct: bool
    is_metal: bool
    efermi: float
    
    def __str__(self) -> str:
        if self.is_metal:
            return f"Metallic (Fermi energy: {self.efermi:.4f} eV)"
        
        lines = [
            f"Band gap: {self.band_gap:.4f} eV ({'direct' if self.is_direct else 'indirect'})",
            f"VBM: {self.vbm_info.energy:.4f} eV at k-point {self.vbm_info.kpoint_coords}",
            f"CBM: {self.cbm_info.energy:.4f} eV at k-point {self.cbm_info.kpoint_coords}",
        ]
        return "\n".join(lines)


@dataclass
class EffectiveMassResult:
    """Effective mass calculation result.
    
    Attributes:
        effective_mass: The effective mass object with detailed data
        temperature: Temperature used in calculation (K)
        concentrations: Carrier concentrations used
    """
    effective_mass: any
    temperature: float
    concentrations: List[float]
    
    def save_json(self, filename: Union[str, Path] = "effective_mass.json") -> None:
        """Save effective mass data as JSON file."""
        self.effective_mass.to_json_file(str(filename))
    
    def __str__(self) -> str:
        return str(self.effective_mass)


def get_band_edge_properties(
    vasprun: Union[str, Path],
    outcar: Union[str, Path] = "OUTCAR",
    integer_criterion: float = 0.1
) -> BandEdgeResult:
    """
    Get band edge properties from VASP output.
    
    Args:
        vasprun: Path to vasprun.xml file
        outcar: Path to OUTCAR file
        integer_criterion: Criterion for integer occupation
    
    Returns:
        BandEdgeResult with band edge information
    
    Example:
        >>> from vise.api.band_edge import get_band_edge_properties
        >>>
        >>> result = get_band_edge_properties("vasprun.xml", "OUTCAR")
        >>> if not result.is_metal:
        ...     print(f"Band gap: {result.band_gap:.2f} eV")
        ...     print(f"VBM at: {result.vbm_info.kpoint_coords}")
    """
    from pymatgen.io.vasp import Vasprun, Outcar
    
    vasprun_obj = Vasprun(str(vasprun))
    outcar_obj = Outcar(str(outcar))
    
    band_edge = VaspBandEdgeProperties(
        vasprun_obj, outcar_obj, 
        integer_criterion=integer_criterion
    )
    
    if band_edge.band_gap is None:
        return BandEdgeResult(
            band_gap=None,
            vbm_info=None,
            cbm_info=None,
            is_direct=False,
            is_metal=True,
            efermi=vasprun_obj.efermi
        )
    
    vbm = band_edge.vbm_info
    cbm = band_edge.cbm_info
    
    vbm_info = KpointInfo(
        kpoint_index=vbm.kpoint_index,
        kpoint_coords=tuple(vbm.kpoint_coords),
        energy=vbm.energy
    )
    
    cbm_info = KpointInfo(
        kpoint_index=cbm.kpoint_index,
        kpoint_coords=tuple(cbm.kpoint_coords),
        energy=cbm.energy
    )
    
    is_direct = vbm.kpoint_index == cbm.kpoint_index
    
    return BandEdgeResult(
        band_gap=band_edge.band_gap,
        vbm_info=vbm_info,
        cbm_info=cbm_info,
        is_direct=is_direct,
        is_metal=False,
        efermi=vasprun_obj.efermi
    )


def calculate_effective_mass(
    vasprun: Union[str, Path],
    outcar: Union[str, Path] = "OUTCAR",
    temperature: float = 300.0,
    concentrations: Optional[List[float]] = None
) -> EffectiveMassResult:
    """
    Calculate effective mass from VASP output.
    
    Requires BoltzTrap2 to be installed.
    
    Args:
        vasprun: Path to vasprun.xml file
        outcar: Path to OUTCAR file
        temperature: Temperature in Kelvin
        concentrations: List of carrier concentrations
    
    Returns:
        EffectiveMassResult with effective mass data
    
    Example:
        >>> from vise.api.band_edge import calculate_effective_mass
        >>>
        >>> result = calculate_effective_mass(
        ...     "vasprun.xml", "OUTCAR",
        ...     temperature=300, 
        ...     concentrations=[1e18, 1e19, 1e20]
        ... )
        >>> print(result)
        >>> result.save_json()
    
    Raises:
        ImportError: If BoltzTrap2 is not installed
        ValueError: If material is metallic
    """
    from pymatgen.io.vasp import Vasprun, Outcar
    
    if concentrations is None:
        concentrations = [1e18, 1e19, 1e20]
    
    vasprun_obj = Vasprun(str(vasprun))
    outcar_obj = Outcar(str(outcar))
    
    band_edge_prop = VaspBandEdgeProperties(vasprun_obj, outcar_obj)
    
    try:
        vbm, cbm = band_edge_prop.vbm_cbm
    except TypeError:
        raise ValueError(
            "Band gap does not exist, so not suited for effective "
            "mass calculation."
        )
    
    effective_mass = make_effective_mass(
        vasprun_obj,
        temperature,
        concentrations,
        vbm, cbm
    )
    
    return EffectiveMassResult(
        effective_mass=effective_mass,
        temperature=temperature,
        concentrations=concentrations
    )
