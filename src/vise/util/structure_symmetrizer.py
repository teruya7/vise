# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Crystal structure symmetry analysis utilities.

This module provides the StructureSymmetrizer class for analyzing crystal
structure symmetry using spglib and seekpath libraries. It handles primitive
cell detection, space group identification, and k-point path generation.
"""

from collections import defaultdict
from dataclasses import dataclass
from itertools import groupby
from typing import Dict, List, Optional, Tuple

import numpy as np
import seekpath
import spglib
from monty.json import MSONable
from more_itertools import consecutive_groups
from pymatgen.core import Element, Structure
from spglib import SpglibDataset
from tabulate import tabulate

from vise.defaults import defaults
from vise.error import ViseError
from vise.util.bravais_lattice import BravaisLattice
from vise.util.centering import Centering
from vise.util.logger import get_logger

logger = get_logger(__name__)

# Type alias for spglib cell format
Cell = Tuple[List[List[float]], List[List[float]], List[int]]

# Mapping of point group symbols to number of symmetry operations
POINT_GROUP_ORDER: Dict[str, int] = {
    "1": 1, "-1": 2, "2": 2, "m": 2, "2/m": 4,
    "222": 4, "2mm": 4, "m2m": 4, "mm2": 4, "mmm": 8,
    "4": 4, "-4": 4, "4/m": 8, "422": 8, "4mm": 8,
    "-4m2": 8, "-42m": 8, "4/mmm": 16,
    "3": 3, "-3": 6, "32": 6, "3m": 6, "-3m": 12,
    "6": 6, "-6": 6, "6/m": 12, "622": 12, "6mm": 12,
    "-6m2": 12, "6/mmm": 24,
    "23": 12, "m3": 24, "m-3": 24, "432": 24, "-43m": 24, "m-3m": 48,
}


def cell_to_structure(cell: Cell) -> Structure:
    """Convert spglib cell format to pymatgen Structure.

    Args:
        cell: Tuple of (lattice_matrix, frac_coords, atomic_numbers).

    Returns:
        pymatgen Structure object.
    """
    species = [Element.from_Z(z) for z in cell[2]]
    return Structure(lattice=cell[0], coords=cell[1], species=species)


def structure_to_cell(structure: Structure) -> Cell:
    """Convert pymatgen Structure to spglib cell format.

    Args:
        structure: pymatgen Structure object.

    Returns:
        Tuple of (lattice_matrix, frac_coords, atomic_numbers).
    """
    lattice = structure.lattice.matrix.tolist()
    frac_coords = structure.frac_coords.tolist()
    z_numbers = [site.specie.Z for site in structure.sites]
    return (lattice, frac_coords, z_numbers)


class StructureSymmetrizer:
    """Symmetry analyzer for crystal structures.

    Uses spglib for symmetry detection and seekpath for band path generation.
    Provides access to primitive/conventional cells, space group information,
    and irreducible k-points.

    Attributes:
        structure: The input crystal structure.
        symprec: Symmetry tolerance for length (Angstrom).
        angle_tolerance: Symmetry tolerance for angles (degrees).
        time_reversal: Whether to consider time reversal symmetry.
        ref_distance: Reference distance for band path k-point sampling.
    """

    def __init__(
        self,
        structure: Structure,
        symprec: float = defaults.symmetry_length_tolerance,
        angle_tolerance: float = defaults.symmetry_angle_tolerance,
        time_reversal: bool = True,
        band_mesh_distance: float = defaults.band_mesh_distance,
    ) -> None:
        """Initialize the symmetry analyzer.

        Args:
            structure: Crystal structure to analyze.
            symprec: Length tolerance for symmetry detection.
            angle_tolerance: Angle tolerance for symmetry detection.
            time_reversal: Consider time reversal symmetry.
            band_mesh_distance: K-point spacing for band paths.
        """
        self.structure = structure.copy()

        if structure.site_properties:
            logger.warning(
                f"Site properties {list(structure.site_properties.keys())} "
                f"will be removed in primitive and conventional structures."
            )

        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        self.time_reversal = time_reversal
        self.ref_distance = band_mesh_distance

        # Build spglib cell
        lattice_matrix = structure.lattice.matrix
        positions = structure.frac_coords.tolist()
        atomic_numbers = [site.specie.number for site in structure.sites]
        self.cell = (lattice_matrix, positions, atomic_numbers)

        # Lazily evaluated properties
        self._spglib_sym_data: Optional[SpglibDataset] = None
        self._primitive: Optional[Structure] = None
        self._second_primitive: Optional[Structure] = None
        self._seekpath_data: Optional[dict] = None

    def __repr__(self) -> str:
        """Return a formatted string representation."""
        sym_data = self.spglib_sym_data
        is_primitive = self.structure == self.primitive

        lines = [
            f"Symprec: {self.symprec}",
            f"Angle tolerance: {self.angle_tolerance}",
            f"Space group: {sym_data['international']}",
            f"Is primitive: {is_primitive}",
        ]

        # Build site information table
        site_headers = ["site", "wyckoff", "site sym", "equiv sites"]
        site_rows = [
            [
                name,
                site.wyckoff_letter,
                site.site_symmetry,
                site.pprint_equiv_atoms,
            ]
            for name, site in self.sites.items()
        ]
        lines.append(tabulate(site_rows, headers=site_headers))

        return "\n".join(lines)

    @property
    def spglib_sym_data(self) -> SpglibDataset:
        """Get spglib symmetry dataset (lazily evaluated)."""
        if self._spglib_sym_data is None:
            self._spglib_sym_data = spglib.get_symmetry_dataset(
                self.cell, self.symprec, self.angle_tolerance
            )
            if self._spglib_sym_data is None:
                raise ValueError("Spglib could not determine symmetry dataset")
        return self._spglib_sym_data

    @property
    def conventional(self) -> Structure:
        """Get the conventional unit cell."""
        center = Centering.from_string(self.centering)
        return self.primitive * center.primitive_to_conv

    @property
    def primitive(self) -> Structure:
        """Get the primitive unit cell (lazily evaluated)."""
        if self._primitive is not None:
            return self._primitive

        primitive_cell = spglib.find_primitive(
            self.cell, symprec=self.symprec, angle_tolerance=self.angle_tolerance
        )
        if primitive_cell is None:
            raise ViseSymmetryError(
                "Spglib couldn't find the primitive cell. "
                "Try adjusting symprec and/or angle_tolerance."
            )

        # Handle spglib cyclic behavior by running twice
        second_primitive = cell_to_structure(
            spglib.find_primitive(
                primitive_cell,
                symprec=self.symprec,
                angle_tolerance=self.angle_tolerance,
            )
        ).get_sorted_structure()

        primitive = cell_to_structure(primitive_cell)

        if primitive != second_primitive:
            if _first_structure_is_primitive(primitive, second_primitive):
                self._primitive = primitive
                self._second_primitive = second_primitive
            else:
                self._primitive = second_primitive
                self._second_primitive = primitive
        else:
            self._primitive = primitive

        return self._primitive

    @property
    def second_primitive(self) -> Optional[Structure]:
        """Get alternate primitive cell if spglib behavior was cyclic."""
        return self._second_primitive

    def find_seekpath_data(self) -> None:
        """Calculate seekpath band path data."""
        logger.info(
            f"Band mesh distance is set to {self.ref_distance}. "
            f"Use band_ref_dist option to change it."
        )
        cell = structure_to_cell(self.primitive)
        self._seekpath_data = seekpath.get_explicit_k_path_orig_cell(
            structure=cell,
            symprec=self.symprec,
            angle_tolerance=self.angle_tolerance,
            with_time_reversal=self.time_reversal,
            reference_distance=self.ref_distance,
        )

    @property
    def sg_number(self) -> int:
        """Space group number."""
        return self.spglib_sym_data.number

    @property
    def space_group(self) -> str:
        """International space group symbol."""
        return self.spglib_sym_data.international

    @property
    def point_group(self) -> str:
        """Point group symbol."""
        return self.spglib_sym_data.pointgroup

    @property
    def seekpath_data(self) -> dict:
        """Seekpath band path data (lazily evaluated)."""
        if self._seekpath_data is None:
            self.find_seekpath_data()
        return self._seekpath_data

    @property
    def is_primitive_lattice_changed(self) -> bool:
        """Check if input structure differs from primitive cell."""
        return self.structure.lattice != self.primitive.lattice

    def irreducible_kpoints(
        self,
        num_kpt_list: List[int],
        kpt_shift: List[float],
    ) -> List[Tuple[List[float], int]]:
        """Calculate irreducible k-points for the given mesh.

        Args:
            num_kpt_list: Number of k-points along each axis [nx, ny, nz].
            kpt_shift: K-point mesh shift in VASP convention [0-0.5].

        Returns:
            List of (kpoint, weight) tuples where kpoint is in
            fractional coordinates.
        """
        mapping, integer_grid_points = spglib.get_ir_reciprocal_mesh(
            mesh=num_kpt_list,
            cell=self.cell,
            is_shift=np.array(kpt_shift) * 2,
            is_time_reversal=self.time_reversal,
            symprec=self.symprec,
        )

        results = []
        shift = np.array(kpt_shift)
        for repr_index, count in zip(*np.unique(mapping, return_counts=True)):
            grid_point = integer_grid_points[repr_index] + shift
            kpoint = (grid_point / num_kpt_list).tolist()
            results.append((kpoint, count))

        return results

    def grouped_atom_indices(self) -> Dict[str, List[int]]:
        """Get atom indices grouped by element and Wyckoff position."""
        return {
            f"{name}_{site.wyckoff_letter}": site.equivalent_atoms
            for name, site in self.sites.items()
        }

    @property
    def sites(self) -> Dict[str, "Site"]:
        """Get site information grouped by symmetry equivalence."""
        wyckoffs = self.spglib_sym_data.wyckoffs
        equivalent_atoms = self.spglib_sym_data.equivalent_atoms
        site_symmetries = self.spglib_sym_data.site_symmetry_symbols

        equiv_indices = sorted(
            enumerate(equivalent_atoms), key=lambda x: x[1]
        )

        result: Dict[str, Site] = {}
        element_counts: Dict[str, int] = defaultdict(int)

        for _, equiv_sites in groupby(equiv_indices, lambda x: x[1]):
            equiv_site_list = list(equiv_sites)
            repr_idx = equiv_site_list[0][0]
            element = self.structure[repr_idx].specie.name

            element_counts[element] += 1
            name = f"{element}{element_counts[element]}"

            result[name] = Site(
                element=element,
                wyckoff_letter=wyckoffs[repr_idx],
                site_symmetry=site_symmetries[repr_idx],
                equivalent_atoms=[s[0] for s in equiv_site_list],
            )

        return result

    @property
    def bravais(self) -> BravaisLattice:
        """Get Bravais lattice type."""
        return BravaisLattice.from_sg_num(self.spglib_sym_data.number)

    @property
    def centering(self) -> str:
        """Get lattice centering type letter."""
        return self.spglib_sym_data.international[0]


def _first_structure_is_primitive(
    structure1: Structure, structure2: Structure
) -> bool:
    """Determine which of two structures is the 'primary' primitive.

    Compares fractional coordinates lexicographically to pick a consistent
    choice when spglib cycles between equivalent primitive cells.
    """
    for s1, s2 in zip(structure1, structure2):
        for c1, c2 in zip(s1.frac_coords, s2.frac_coords):
            if c1 < c2 - 1e-5:
                return True
            elif c2 < c1 - 1e-5:
                return False
    raise ViseSymmetryError("Cannot uniquely determine primitive cell.")


@dataclass(frozen=True)
class Site(MSONable):
    """Crystallographic site information.

    Attributes:
        element: Element symbol.
        wyckoff_letter: Wyckoff position letter.
        site_symmetry: Site symmetry symbol.
        equivalent_atoms: List of atom indices at equivalent positions.
    """

    element: str
    wyckoff_letter: str
    site_symmetry: str
    equivalent_atoms: List[int]

    @property
    def pprint_equiv_atoms(self) -> str:
        """Format equivalent atom indices for display.

        Consecutive ranges are shown as 'start..end'.
        """
        parts = []
        for consecutive in consecutive_groups(self.equivalent_atoms):
            indices = list(consecutive)
            if len(indices) >= 3:
                parts.append(f"{indices[0]}..{indices[-1]}")
            else:
                parts.append(" ".join(str(i) for i in indices))
        return " ".join(parts)


def num_symmetry_operation(point_group: str) -> int:
    """Get the order (number of operations) of a point group.

    Args:
        point_group: Hermann-Mauguin point group symbol.
                    Dots are ignored (e.g., "..6" -> "6").

    Returns:
        Number of symmetry operations in the point group.

    Examples:
        >>> num_symmetry_operation("m-3m")
        48
        >>> num_symmetry_operation("..6")
        6
    """
    # Remove orientation dots from the symbol
    cleaned = point_group.replace(".", "")
    return POINT_GROUP_ORDER[cleaned]


# Keep old name for backward compatibility
num_sym_op = POINT_GROUP_ORDER


class ViseSymmetryError(ViseError):
    """Exception raised for symmetry analysis errors."""

    pass


# Backward compatibility alias
first_structure_is_primitive = _first_structure_is_primitive

