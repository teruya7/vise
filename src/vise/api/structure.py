# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Structure analysis API.

This module provides functions for analyzing crystal structure symmetry
without CLI dependencies.

Example:
    >>> from pymatgen.core import Structure
    >>> from vise.api.structure import get_symmetry_info
    >>>
    >>> struct = Structure.from_file("POSCAR")
    >>> info = get_symmetry_info(struct)
    >>> print(f"Space group: {info.space_group_symbol} ({info.space_group_number})")
    >>> print(f"Crystal system: {info.crystal_system}")
"""

from dataclasses import dataclass
from typing import Optional, Tuple

from pymatgen.core import Structure

from vise.defaults import defaults
from vise.util.structure_symmetrizer import StructureSymmetrizer


@dataclass
class SymmetryInfo:
    """Structure symmetry information.
    
    Attributes:
        space_group_symbol: International space group symbol (e.g., "Fm-3m")
        space_group_number: Space group number (1-230)
        point_group: Point group symbol
        crystal_system: Crystal system name (e.g., "cubic")
        bravais_lattice: Bravais lattice name
        centering: Centering type (P, I, F, C, R, A, B)
        volume: Unit cell volume in Å³
        lattice_abc: Lattice parameters (a, b, c) in Å
        lattice_angles: Lattice angles (alpha, beta, gamma) in degrees
        is_primitive: Whether the input structure is already primitive
    """
    space_group_symbol: str
    space_group_number: int
    point_group: str
    crystal_system: str
    bravais_lattice: str
    centering: str
    volume: float
    lattice_abc: Tuple[float, float, float]
    lattice_angles: Tuple[float, float, float]
    is_primitive: bool


def get_symmetry_info(
    structure: Structure,
    symprec: Optional[float] = None,
    angle_tolerance: Optional[float] = None
) -> SymmetryInfo:
    """
    Analyze structure symmetry and return detailed symmetry information.
    
    Args:
        structure: pymatgen Structure object to analyze
        symprec: Length tolerance in Å for symmetry analysis.
            Default: defaults.symmetry_length_tolerance (0.01)
        angle_tolerance: Angle tolerance in degrees for symmetry analysis.
            Default: defaults.symmetry_angle_tolerance (5.0)
    
    Returns:
        SymmetryInfo dataclass containing space group, point group,
        crystal system, and lattice information
    
    Example:
        >>> from pymatgen.core import Structure
        >>> from vise.api.structure import get_symmetry_info
        >>>
        >>> struct = Structure.from_file("POSCAR")
        >>> info = get_symmetry_info(struct)
        >>> print(info.space_group_symbol)  # e.g., "Fm-3m"
        >>> print(info.space_group_number)  # e.g., 225
        >>> print(info.crystal_system)      # e.g., "cubic"
    """
    symprec = symprec if symprec is not None else defaults.symmetry_length_tolerance
    angle_tolerance = angle_tolerance if angle_tolerance is not None else defaults.symmetry_angle_tolerance
    
    symmetrizer = StructureSymmetrizer(
        structure, 
        symprec=symprec, 
        angle_tolerance=angle_tolerance
    )
    
    sym_data = symmetrizer.spglib_sym_data
    is_primitive = structure == symmetrizer.primitive
    
    # Extract crystal system from BravaisLattice name (first letter)
    bravais_name = symmetrizer.bravais.name
    crystal_system_map = {
        'a': 'triclinic',
        'm': 'monoclinic',
        'o': 'orthorhombic',
        't': 'tetragonal',
        'h': 'hexagonal',  # includes trigonal (hR)
        'c': 'cubic'
    }
    crystal_system = crystal_system_map.get(bravais_name[0], 'unknown')
    
    return SymmetryInfo(
        space_group_symbol=sym_data.international,
        space_group_number=sym_data.number,
        point_group=sym_data.pointgroup,
        crystal_system=crystal_system,
        bravais_lattice=bravais_name,
        centering=symmetrizer.centering,
        volume=structure.volume,
        lattice_abc=structure.lattice.abc,
        lattice_angles=structure.lattice.angles,
        is_primitive=is_primitive
    )


def get_primitive(
    structure: Structure,
    symprec: Optional[float] = None,
    angle_tolerance: Optional[float] = None
) -> Structure:
    """
    Get the primitive cell of a crystal structure.
    
    Args:
        structure: pymatgen Structure object
        symprec: Length tolerance in Å for symmetry analysis.
            Default: defaults.symmetry_length_tolerance (0.01)
        angle_tolerance: Angle tolerance in degrees for symmetry analysis.
            Default: defaults.symmetry_angle_tolerance (5.0)
    
    Returns:
        Primitive cell as a pymatgen Structure object.
        Returns the input structure if it's already primitive.
    
    Example:
        >>> from pymatgen.core import Structure
        >>> from vise.api.structure import get_primitive
        >>>
        >>> conventional = Structure.from_file("POSCAR")
        >>> primitive = get_primitive(conventional)
        >>> print(f"Primitive cell has {len(primitive)} atoms")
    """
    symprec = symprec if symprec is not None else defaults.symmetry_length_tolerance
    angle_tolerance = angle_tolerance if angle_tolerance is not None else defaults.symmetry_angle_tolerance
    
    symmetrizer = StructureSymmetrizer(
        structure, 
        symprec=symprec, 
        angle_tolerance=angle_tolerance
    )
    return symmetrizer.primitive


def get_conventional(
    structure: Structure,
    symprec: Optional[float] = None,
    angle_tolerance: Optional[float] = None
) -> Structure:
    """
    Get the conventional cell of a crystal structure.
    
    Args:
        structure: pymatgen Structure object
        symprec: Length tolerance in Å for symmetry analysis.
            Default: defaults.symmetry_length_tolerance (0.01)  
        angle_tolerance: Angle tolerance in degrees for symmetry analysis.
            Default: defaults.symmetry_angle_tolerance (5.0)
    
    Returns:
        Conventional cell as a pymatgen Structure object.
        Returns the input structure if it's already conventional.
    
    Example:
        >>> from pymatgen.core import Structure
        >>> from vise.api.structure import get_conventional
        >>>
        >>> primitive = Structure.from_file("POSCAR")
        >>> conventional = get_conventional(primitive)
        >>> print(f"Conventional cell has {len(conventional)} atoms")
    """
    symprec = symprec if symprec is not None else defaults.symmetry_length_tolerance
    angle_tolerance = angle_tolerance if angle_tolerance is not None else defaults.symmetry_angle_tolerance
    
    symmetrizer = StructureSymmetrizer(
        structure, 
        symprec=symprec, 
        angle_tolerance=angle_tolerance
    )
    return symmetrizer.conventional
