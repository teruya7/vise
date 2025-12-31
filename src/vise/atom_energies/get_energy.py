# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Utilities for extracting atomic energies from VASP calculations.

This module provides functions to extract converged atomic energies
from VASP output files for building reference energy databases.
"""

import os
import warnings
from pathlib import Path
from xml.etree.ElementTree import ParseError

from pymatgen.core import Element
from pymatgen.io.vasp import Outcar, Vasprun
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

from vise.atom_energies.make_atom_vasp_set import is_target_element
from vise.util.logger import get_logger

logger = get_logger(__name__)

# Suppress POTCAR warnings for cleaner output
warnings.simplefilter("ignore", UnknownPotcarWarning)


def make_energy_yaml() -> None:
    """Extract atomic energies from subdirectories and print YAML format.

    Scans current directory for element-named subdirectories, extracts
    converged energies from VASP calculations, and prints in YAML format.

    Output format:
        H:  -3.12345678
        He: -0.00012345
        ...

    Notes:
        - Skips elements with unconverged calculations
        - Skips elements with parse errors in vasprun.xml
        - Only processes target elements (Z<=58 or 71<=Z<=83)
    """
    dirs = [f.name for f in os.scandir(".") if f.is_dir()]

    for elem in list(Element):
        elem_str = str(elem)

        if elem_str not in dirs or not is_target_element(elem_str):
            continue

        try:
            vasprun = Vasprun(Path(elem_str) / "vasprun.xml")

            if not vasprun.converged_electronic:
                logger.warning(f"Calculation for {elem_str} not converged.")
                continue
        except ParseError:
            logger.warning(f"Failed to parse vasprun.xml for {elem_str}.")
            continue

        outcar = Outcar(Path(elem_str) / "OUTCAR")
        print(f"{elem_str + ':':<3} {outcar.final_energy:11.8f}")


if __name__ == "__main__":
    make_energy_yaml()
