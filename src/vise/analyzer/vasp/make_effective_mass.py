# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Effective mass calculation from VASP using BoltzTraP2.

This module provides utilities for calculating effective masses from
VASP band structure calculations using the BoltzTraP2 transport code.
"""

from copy import deepcopy
from typing import List

import numpy as np
from pymatgen.electronic_structure.boltztrap import BoltztrapError
from pymatgen.io.vasp import Vasprun

from vise.analyzer.effective_mass import EffectiveMass


def make_effective_mass(
    vasprun: Vasprun,
    temp: float,
    concentrations: List[float],
    vbm: float,
    cbm: float,
) -> EffectiveMass:
    """Calculate effective masses using BoltzTraP2.

    Uses BoltzTraP2 to interpolate bands and compute the effective
    mass tensor at specified temperature and carrier concentrations.

    Args:
        vasprun: Parsed vasprun.xml with band structure.
        temp: Temperature in Kelvin.
        concentrations: List of carrier concentrations in cm^-3.
        vbm: VBM energy in eV.
        cbm: CBM energy in eV.

    Returns:
        EffectiveMass with hole and electron mass tensors.

    Raises:
        ImportError: If BoltzTraP2 is not installed.

    Examples:
        >>> em = make_effective_mass(vasprun, 300, [1e18, 1e19], -5.0, -3.5)
        >>> em.average_mass("p", 1e18)
        0.35
    """
    try:
        from pymatgen.electronic_structure.boltztrap2 import (
            BztInterpolator,
            BztTransportProperties,
            VasprunBSLoader,
        )
    except BoltztrapError:
        raise ImportError("Calculating effective mass requires BoltzTraP2")

    # Deep copy to avoid side effects from Fermi level modification
    vasprun = deepcopy(vasprun)

    # Set Fermi level to mid-gap to fix pymatgen band edge assignment
    vasprun.efermi = (cbm + vbm) / 2

    # Load and interpolate band structure
    loader = VasprunBSLoader(vasprun)
    energy_range = (cbm - vbm) / 2 + 2.0
    interpolator = BztInterpolator(loader, energy_range=energy_range)

    # Compute transport properties
    transport = BztTransportProperties(interpolator, temp_r=np.array([temp]))
    transport.compute_properties_doping(concentrations)

    return EffectiveMass(
        p=transport.Effective_mass_doping["p"].tolist()[0],
        n=transport.Effective_mass_doping["n"].tolist()[0],
        temperature=temp,
        concentrations=concentrations,
    )