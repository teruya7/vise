# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Dielectric function extraction from VASP.

This module provides utilities for creating DieleFuncData from VASP
optical calculations (LOPTICS = .TRUE.).
"""

from typing import Optional

import numpy as np
from pymatgen.io.vasp import Outcar, Vasprun

from vise.analyzer.dielectric_function import DieleFuncData, kramers_kronig_trans
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties


def make_diele_func(
    vasprun: Vasprun,
    outcar: Outcar,
    use_vasp_real: bool = True,
    ita: float = 0.01,
) -> DieleFuncData:
    """Create DieleFuncData from VASP optical calculation.

    Args:
        vasprun: Parsed vasprun.xml with dielectric data.
        outcar: Parsed OUTCAR for band edge determination.
        use_vasp_real: If True, use VASP's real part directly.
                      If False, compute via Kramers-Kronig transform.
        ita: Broadening parameter for Kramers-Kronig (only if use_vasp_real=False).

    Returns:
        DieleFuncData with dielectric function components.
    """
    energies, real, imag = vasprun.dielectric_data["density"]
    imag = np.array(imag)

    if use_vasp_real:
        real = np.array(real)
        # Transpose to get [component][energy] ordering
        real, imag = real.T, imag.T
    else:
        # When CSHIFT = 0.0, first components can be unphysical (99999.0)
        # See vasprun.xml structure:
        # <dielectricfunction>
        #   <imag>
        #     <field>xx</field>...
        #     <r> 0.0 0.0 0.0 -0.0 99999.0 99999.0 99999.0 </r>
        imag[0] = 0.0
        imag = imag.T
        real = kramers_kronig_trans(imag, energies, ita)

    # Add isotropic average
    real = np.vstack([real, _make_average(real)])
    imag = np.vstack([imag, _make_average(imag)])

    # Get band gap from VBM/CBM
    band_gap = VaspBandEdgeProperties(vasprun, outcar).band_gap

    return DieleFuncData(
        energies=energies,
        directions=["xx", "yy", "zz", "xy", "yz", "xz", "ave"],
        diele_func_real=real.tolist(),
        diele_func_imag=imag.tolist(),
        band_gap=band_gap,
    )


def _make_average(vals: np.ndarray) -> np.ndarray:
    """Calculate isotropic average from diagonal components.

    Args:
        vals: Array with shape [6, n_energy] (xx, yy, zz, xy, yz, xz).

    Returns:
        Average of diagonal components (xx, yy, zz).
    """
    return np.average(vals[:3, :], axis=0)


# Backward compatibility alias
make_average = _make_average
