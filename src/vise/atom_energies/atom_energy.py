# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
"""Atomic reference energies.

This module provides atomic reference energies from Materials Project
calculations, used for computing formation energies.
"""

from pathlib import Path

from monty.serialization import loadfn

from vise.util.logger import get_logger

logger = get_logger(__name__)

# Path to data files
_DATA_DIR = Path(__file__).parent

# Materials Project atomic energies (eV/atom)
mp_energies = loadfn(_DATA_DIR / "mp_atom_energy.yaml")
