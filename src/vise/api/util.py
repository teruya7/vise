# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Utility functions API.

This module provides utility functions for various VASP-related tasks:
- Creating atom POSCARs for reference calculations
- Handling volumetric data (CHGCAR, etc.)
- Phonon calculation setup

Example:
    >>> from vise.api.util import make_atom_poscars
    >>>
    >>> make_atom_poscars(["Si", "O"], output_dir="./atoms")
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

from pymatgen.core import Element, Structure


def make_atom_poscars(
    elements: Optional[List[Union[str, Element]]] = None,
    output_dir: Union[str, Path] = "."
) -> List[Path]:
    """
    Create POSCAR files for isolated atom calculations.
    
    Args:
        elements: List of elements. If None, all elements are created.
        output_dir: Directory where atom calculation directories are created.
    
    Returns:
        List of created directory paths.
    
    Example:
        >>> from vise.api.util import make_atom_poscars
        >>>
        >>> # Create poscars for specific elements
        >>> dirs = make_atom_poscars(["Si", "O"], output_dir="./atoms")
        >>> print(f"Created {len(dirs)} directories")
        >>>
        >>> # Create for all elements
        >>> dirs = make_atom_poscars()
    """
    from vise.atom_energies.make_atom_vasp_set import make_atom_poscar_dirs

    output_path = Path(output_dir)

    # Convert strings to Element objects if needed
    if elements:
        element_objs = [
            Element(e) if isinstance(e, str) else e
            for e in elements
        ]
    else:
        element_objs = None

    make_atom_poscar_dirs(output_path, element_objs)

    # Return list of created directories
    if element_objs:
        return [output_path / str(e) for e in element_objs]
    else:
        # All elements
        return list(output_path.iterdir())


def make_spin_decomposed_volumetric_files(
    chgcar_file: Union[str, Path],
    output_prefix: Optional[str] = None
) -> tuple:
    """
    Create spin-decomposed volumetric files from CHGCAR.
    
    Args:
        chgcar_file: Path to CHGCAR-type file with spin data.
        output_prefix: Prefix for output files. Default: same as input.
    
    Returns:
        Tuple of (spin_up_file, spin_down_file) paths.
    
    Example:
        >>> from vise.api.util import make_spin_decomposed_volumetric_files
        >>>
        >>> up, down = make_spin_decomposed_volumetric_files("CHGCAR")
        >>> print(f"Created: {up}, {down}")
    """
    from pymatgen.io.vasp import Chgcar

    from vise.analyzer.vasp.handle_volumetric_data import make_spin_charges

    chgcar = Chgcar.from_file(str(chgcar_file))

    prefix = output_prefix or str(chgcar_file)
    up_file = f"{prefix}_up"
    down_file = f"{prefix}_down"

    for c, (spin, filename) in zip(
        make_spin_charges(chgcar),
        [("up", up_file), ("down", down_file)]
    ):
        c.write_file(filename)

    return Path(up_file), Path(down_file)


def make_light_weight_volumetric_data(
    volumetric_file: Union[str, Path],
    output_filename: Optional[Union[str, Path]] = None,
    border_fractions: Optional[List[float]] = None
) -> Path:
    """
    Create a lightweight version of volumetric data file.
    
    Args:
        volumetric_file: Path to CHGCAR-type volumetric file.
        output_filename: Output filename. Default: {input}_lw
        border_fractions: Fractions for isosurfaces.
    
    Returns:
        Path to created lightweight file.
    
    Example:
        >>> from vise.api.util import make_light_weight_volumetric_data
        >>>
        >>> lw_file = make_light_weight_volumetric_data("CHGCAR")
        >>> print(f"Created: {lw_file}")
    """
    from pymatgen.io.vasp import Chgcar

    from vise.analyzer.vasp.handle_volumetric_data import (
        default_border_fractions,
        light_weight_vol_text,
    )

    chgcar = Chgcar.from_file(str(volumetric_file))

    output_path = Path(output_filename or f"{volumetric_file}_lw")
    fractions = border_fractions or default_border_fractions

    lw_text = light_weight_vol_text(chgcar, fractions)
    output_path.write_text(lw_text)

    return output_path


def create_vesta_file(
    volumetric_file: Union[str, Path],
    output_vesta_file: Union[str, Path],
    border_fractions: Optional[List[float]] = None,
    original_vesta_file: Optional[Union[str, Path]] = None,
    lw_volumetric_file: Optional[Union[str, Path]] = None
) -> Path:
    """
    Create a VESTA file with volumetric data and isosurfaces.
    
    Args:
        volumetric_file: Path to CHGCAR-type volumetric file.
        output_vesta_file: Output VESTA filename.
        border_fractions: Fractions for isosurfaces.
        original_vesta_file: Original VESTA file to use as base.
        lw_volumetric_file: Lightweight volumetric file path to reference.
    
    Returns:
        Path to created VESTA file.
    
    Example:
        >>> from vise.api.util import create_vesta_file
        >>>
        >>> vesta = create_vesta_file("CHGCAR", "structure.vesta")
    """
    from pymatgen.io.vasp import Chgcar

    from vise.analyzer.vasp.handle_volumetric_data import default_border_fractions
    from vise.analyzer.vesta.vesta_file import VestaFile, add_density, calc_isurfs

    chgcar = Chgcar.from_file(str(volumetric_file))
    fractions = border_fractions or default_border_fractions

    is_chg = "CHG" in str(volumetric_file)
    volume = chgcar.structure.volume
    isurfs = calc_isurfs(fractions, is_chg, volume)

    if original_vesta_file:
        vesta_text = Path(original_vesta_file).read_text()
    else:
        vesta_text = VestaFile(chgcar.structure).__repr__()

    vol_filename = str(lw_volumetric_file or volumetric_file)
    to_vesta_text = add_density(vesta_text, isurfs, vol_filename)

    output_path = Path(output_vesta_file)
    output_path.write_text(to_vesta_text)

    return output_path


@dataclass
class PhonopySetup:
    """Phonopy calculation setup result.
    
    Attributes:
        supercell: Supercell structure
        phonopy_input_file: Path to phonopy_input.json
    """
    supercell: Structure
    phonopy_input_file: Path


def create_phonon_setup(
    unitcell: Union[str, Path, Structure],
    supercell_matrix: List[int],
    output_dir: Union[str, Path] = "."
) -> PhonopySetup:
    """
    Create phonon calculation setup.
    
    Args:
        unitcell: Path to POSCAR or Structure object.
        supercell_matrix: Supercell matrix (3 integers for diagonal,
            or 9 integers for full matrix).
        output_dir: Output directory.
    
    Returns:
        PhonopySetup with supercell and input file path.
    
    Example:
        >>> from vise.api.util import create_phonon_setup
        >>>
        >>> setup = create_phonon_setup("POSCAR", [2, 2, 2])
        >>> print(f"Supercell atoms: {len(setup.supercell)}")
    """
    from vise.util.phonopy.phonopy_input import make_phonopy_input
    from vise.util.structure_handler import sanitize_matrix

    if isinstance(unitcell, (str, Path)):
        unitcell = Structure.from_file(str(unitcell))

    output_path = Path(output_dir)

    phonopy_input = make_phonopy_input(
        unitcell=unitcell,
        supercell_matrix=sanitize_matrix(supercell_matrix)
    )

    # Save supercell POSCAR
    poscar_path = output_path / "POSCAR"
    phonopy_input.supercell.to(filename=str(poscar_path))

    # Save phonopy input
    json_path = output_path / "phonopy_input.json"
    phonopy_input.to_json_file(str(json_path))

    return PhonopySetup(
        supercell=phonopy_input.supercell,
        phonopy_input_file=json_path
    )


def analyze_phonon(
    phonopy_input_file: Union[str, Path],
    vasprun_file: Union[str, Path],
    output_filename: str = "phonon_band.pdf"
) -> Path:
    """
    Analyze phonon calculation and create band plot.
    
    Args:
        phonopy_input_file: Path to phonopy_input.json.
        vasprun_file: Path to vasprun.xml from force calculation.
        output_filename: Output PDF filename.
    
    Returns:
        Path to created plot file.
    
    Example:
        >>> from vise.api.util import analyze_phonon
        >>>
        >>> plot_file = analyze_phonon(
        ...     "phonopy_input.json", 
        ...     "vasprun.xml"
        ... )
    """
    from monty.serialization import loadfn

    phonopy_input = loadfn(str(phonopy_input_file))
    phonopy_input.set_force_constants_from_vasprun(str(vasprun_file))

    phonopy = phonopy_input.to_phonopy
    plt = phonopy.auto_band_structure(plot=True)
    plt.savefig(output_filename)
    plt.close()

    return Path(output_filename)
