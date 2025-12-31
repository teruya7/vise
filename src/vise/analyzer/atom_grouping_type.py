# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Atom grouping utilities for DOS and other analyses.

This module provides different strategies for grouping atoms in crystal
structures, such as by chemical element, by individual atom, or by
symmetrically equivalent sites.
"""

from typing import Dict, List, Optional

from pymatgen.core.structure import Structure

from vise.util.enum import ExtendedEnum
from vise.util.structure_symmetrizer import StructureSymmetrizer


class AtomGroupingType(ExtendedEnum):
    """Atom grouping strategies for analysis.

    Attributes:
        non_equiv_sites: Group by symmetrically non-equivalent sites.
        elements: Group by chemical element type.
        atoms: Group by user-specified atom indices.
    """

    non_equiv_sites = "non_equiv_sites"
    elements = "elements"
    atoms = "atoms"

    def grouped_atom_indices(
        self,
        structure: Structure,
        target: Optional[List[str]] = None,
    ) -> Dict[str, List[int]]:
        """Get atom indices grouped according to this strategy.

        Args:
            structure: Crystal structure to analyze.
            target: For 'atoms' mode, list of comma-separated index strings.
                   For 'elements' mode, list of element symbols.

        Returns:
            Dictionary mapping group names to lists of atom indices.
        """
        if self is self.atoms:
            return group_by_atoms(structure, target)
        elif self is self.elements:
            return group_by_elements(structure, target)
        else:
            return group_by_non_equiv_sites(structure)


def group_by_atoms(
    structure: Structure,
    target: List[str],
) -> Dict[str, List[int]]:
    """Group atoms by user-specified index lists.

    Args:
        structure: Crystal structure (for validation).
        target: List of comma-separated index strings.
               Example: ["0,1,2", "3,4,5"] groups atoms 0-2 and 3-5.

    Returns:
        Dictionary mapping index strings to atom index lists.

    Raises:
        ValueError: If any index exceeds structure size.

    Examples:
        >>> group_by_atoms(structure, ["0,1", "2,3"])
        {'0,1': [0, 1], '2,3': [2, 3]}
    """
    result: Dict[str, List[int]] = {}

    for index_str in target:
        indices = [int(x) for x in index_str.split(",")]
        result[index_str] = indices

    max_index = max(max(indices) for indices in result.values())
    if max_index >= len(structure):
        raise ValueError(
            f"Atom index {max_index} out of range for "
            f"structure with {len(structure)} atoms."
        )

    return result


def group_by_elements(
    structure: Structure,
    target: Optional[List[str]] = None,
) -> Dict[str, List[int]]:
    """Group atoms by chemical element.

    Args:
        structure: Crystal structure to analyze.
        target: Specific elements to include (None for all elements).

    Returns:
        Dictionary mapping element symbols to atom index lists.

    Examples:
        >>> group_by_elements(structure)  # NaCl
        {'Na': [0], 'Cl': [1]}
    """
    elements = target or [str(e) for e in structure.composition.elements]
    result: Dict[str, List[int]] = {}

    for elem in elements:
        result[elem] = [
            i for i, site in enumerate(structure)
            if elem in site
        ]


    return result


def group_by_non_equiv_sites(structure: Structure) -> Dict[str, List[int]]:
    """Group atoms by symmetrically non-equivalent sites.

    Uses spglib symmetry analysis to identify equivalent atoms.

    Args:
        structure: Crystal structure to analyze.

    Returns:
        Dictionary mapping site labels to atom index lists.
    """
    symmetrizer = StructureSymmetrizer(structure)
    return symmetrizer.grouped_atom_indices()
