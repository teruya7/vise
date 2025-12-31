# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Materials Project API.

This module provides functions for interacting with Materials Project
to retrieve structure data.

Example:
    >>> from vise.api.materials_project import get_structure_by_id
    >>>
    >>> structure = get_structure_by_id("mp-149")
    >>> structure.to_file("POSCAR")
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

from pymatgen.core import Structure


@dataclass
class MaterialEntry:
    """Materials Project entry information.
    
    Attributes:
        material_id: Materials Project ID (e.g., "mp-149")
        energy_above_hull: Energy above hull in eV/atom
        space_group: Space group symbol
        band_gap: Band gap in eV
    """
    material_id: str
    energy_above_hull: float
    space_group: str
    band_gap: float


def search_materials(formula: str) -> List[MaterialEntry]:
    """
    Search materials by formula in Materials Project.
    
    Requires MP API key to be configured in ~/.pmgrc.yaml
    
    Args:
        formula: Chemical formula (e.g., "Si", "MgO", "LiFePO4")
    
    Returns:
        List of MaterialEntry sorted by energy above hull
    
    Example:
        >>> from vise.api.materials_project import search_materials
        >>>
        >>> entries = search_materials("TiO2")
        >>> for entry in entries[:5]:
        ...     print(f"{entry.material_id}: {entry.space_group}")
    """
    from mp_api.client import MPRester

    with MPRester() as m:
        candidates = m.summary.search(
            formula=formula,
            fields=["material_id", "energy_above_hull", "symmetry", "band_gap"]
        )

    sorted_candidates = sorted(candidates, key=lambda x: x.energy_above_hull)

    entries = []
    for c in sorted_candidates:
        entries.append(MaterialEntry(
            material_id=str(c.material_id),
            energy_above_hull=c.energy_above_hull,
            space_group=c.symmetry.symbol,
            band_gap=c.band_gap
        ))

    return entries


def get_structure_by_id(
    material_id: str,
    save_poscar: bool = False,
    poscar_filename: Union[str, Path] = "POSCAR"
) -> Structure:
    """
    Get structure from Materials Project by ID.
    
    Args:
        material_id: Materials Project ID (e.g., "mp-149")
        save_poscar: Whether to save structure as POSCAR
        poscar_filename: Output POSCAR filename
    
    Returns:
        pymatgen Structure object
    
    Example:
        >>> from vise.api.materials_project import get_structure_by_id
        >>>
        >>> # Get Si structure
        >>> structure = get_structure_by_id("mp-149")
        >>> print(f"Formula: {structure.composition.reduced_formula}")
        >>>
        >>> # Get and save as POSCAR
        >>> structure = get_structure_by_id("mp-149", save_poscar=True)
    """
    from mp_api.client import MPRester

    structure = MPRester().get_structure_by_material_id(material_id)

    if save_poscar:
        structure.to(fmt="poscar", filename=str(poscar_filename))

    return structure


@dataclass
class MaterialInfo:
    """Detailed material information from Materials Project.
    
    Attributes:
        structure: pymatgen Structure object
        material_id: Materials Project ID
        total_magnetization: Total magnetization
        band_gap: Band gap in eV
    """
    structure: Structure
    material_id: str
    total_magnetization: float
    band_gap: float


def get_material_info(material_id: str) -> MaterialInfo:
    """
    Get detailed material information from Materials Project.
    
    Args:
        material_id: Materials Project ID
    
    Returns:
        MaterialInfo with structure and properties
    
    Example:
        >>> from vise.api.materials_project import get_material_info
        >>>
        >>> info = get_material_info("mp-149")
        >>> print(f"Band gap: {info.band_gap} eV")
        >>> print(f"Magnetization: {info.total_magnetization}")
    """
    from mp_api.client import MPRester

    structure = MPRester().get_structure_by_material_id(material_id)

    query = MPRester().summary.search(
        material_ids=[material_id],
        fields=["total_magnetization", "band_gap"]
    )
    data = query[0]

    return MaterialInfo(
        structure=structure,
        material_id=material_id,
        total_magnetization=data.total_magnetization or 0.0,
        band_gap=data.band_gap
    )


def get_most_stable_structure(
    formula: str,
    save_poscar: bool = False,
    poscar_filename: Union[str, Path] = "POSCAR"
) -> Optional[Structure]:
    """
    Get the most stable structure for a formula from Materials Project.
    
    Args:
        formula: Chemical formula
        save_poscar: Whether to save as POSCAR
        poscar_filename: Output POSCAR filename
    
    Returns:
        Structure of most stable polymorph, or None if not found
    
    Example:
        >>> from vise.api.materials_project import get_most_stable_structure
        >>>
        >>> structure = get_most_stable_structure("TiO2", save_poscar=True)
        >>> if structure:
        ...     print(f"Got {structure.composition.reduced_formula}")
    """
    entries = search_materials(formula)

    if not entries:
        return None

    most_stable_id = entries[0].material_id
    return get_structure_by_id(
        most_stable_id,
        save_poscar=save_poscar,
        poscar_filename=poscar_filename
    )
