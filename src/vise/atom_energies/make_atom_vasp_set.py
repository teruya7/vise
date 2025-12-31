# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""VASP input generation for isolated atom calculations.

This module provides functions to create VASP input files for
calculating atomic reference energies in a large cubic cell.
"""

import warnings
from pathlib import Path
from typing import List, Optional, Union

import yaml
from monty.serialization import loadfn
from pymatgen.core import Element, Lattice, Structure
from pymatgen.io.vasp import Kpoints
from pymatgen.io.vasp.sets import BadInputSetWarning, MPRelaxSet

from vise.input_set.datasets.potcar_set import PotcarSet

# Load atomic magnetization data for spin-polarized calculations
_DATA_DIR = Path(__file__).parent
MAGNETIZATIONS = loadfn(_DATA_DIR / "mp_atom_mag.yaml")

# Suppress MPRelaxSet warnings about isolated atoms
warnings.simplefilter("ignore", BadInputSetWarning)


def is_target_element(elem: Union[Element, str]) -> bool:
    """Check if element should be included in atomic energy calculations.

    Args:
        elem: Element or element symbol.

    Returns:
        True for elements with Z<=58 or 71<=Z<=83 (excludes most lanthanides).
    """
    if isinstance(elem, str):
        elem = Element(elem)
    return elem.Z <= 58 or 71 <= elem.Z <= 83


def make_atom_poscar_dirs(
    path: Path,
    elems: Optional[List[Element]] = None,
) -> None:
    """Create directories with POSCAR files for atomic calculations.

    Creates one directory per element with a POSCAR containing a single
    atom in a 10 Å cubic cell, plus prior_info.yaml for spin settings.

    Args:
        path: Parent directory for element subdirectories.
        elems: Specific elements to include (default: all in magnetization data).
    """
    elements = [str(e) for e in elems] if elems else list(MAGNETIZATIONS.keys())

    for element in elements:
        dir_path = path / element
        dir_path.mkdir()

        # Single atom at cell center
        structure = Structure(
            Lattice.cubic(10),
            [element],
            [[0.5, 0.5, 0.5]],
        )
        structure.to(fmt="POSCAR", filename=str(dir_path / "POSCAR"))

        # Spin-polarized settings with Hund's rule magnetization
        nupdown = MAGNETIZATIONS[element]
        prior_info = {
            "incar": {
                "ISPIN": 2,
                "NUPDOWN": float(nupdown),
                "NELM": 300,
            },
            "is_cluster": True,
        }
        (dir_path / "prior_info.yaml").write_text(yaml.dump(prior_info))


def make_atom_mp_relax_set() -> None:
    """Create Materials Project-style VASP inputs for atomic calculations.

    Creates directories with MPRelaxSet inputs for each element,
    configured for isolated atom calculations (single k-point, etc.).
    """
    for element, potcar in PotcarSet.mp_relax_set.potcar_dict().items():
        if potcar is None:
            continue

        if not is_target_element(element):
            continue

        # Single atom at cell center
        structure = Structure(
            Lattice.cubic(10),
            coords=[[0.5, 0.5, 0.5]],
            species=[element],
        )

        mp_set = MPRelaxSet(
            structure,
            user_kpoints_settings=Kpoints(kpts=((1, 1, 1),)),
            user_incar_settings={
                "ALGO": "D",
                "ISIF": 2,
                "ISMEAR": 0,
                "MAGMOM": {"H": 1.0},
                "NELM": 300,
            },
        )

        Path(element).mkdir()
        mp_set.write_input(element)


if __name__ == "__main__":
    make_atom_mp_relax_set()
