# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Structure and k-point generation for VASP calculations.

This module provides StructureKpointsGenerator for creating properly
configured crystal structures and k-point meshes for various VASP
calculation types.
"""

from math import ceil
from typing import List, Optional, Tuple, Union

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import Kpoints

from vise.defaults import defaults
from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.task import Task
from vise.util.bravais_lattice import BravaisLattice
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


class StructureKpointsGenerator:
    """Generator for structure and k-points based on calculation task.

    Handles symmetry analysis, primitive cell conversion, and k-point
    mesh generation according to the specified calculation task and mode.

    Attributes:
        structure: The (possibly transformed) structure for calculation.
        kpoints: Generated Kpoints object.
        num_kpts: Number of irreducible k-points.
        num_kpt_factor: K-point density multiplier.
        bravais: Bravais lattice type.
    """

    def __init__(
        self,
        initial_structure: Structure,
        task: Task,
        kpt_density: float,
        kpt_mode: Union[KpointsMode, str, None] = None,
        gamma_centered: Optional[bool] = None,
        only_even_num_kpts: bool = False,
        num_kpt_factor: Optional[int] = None,
        band_ref_dist: float = defaults.band_mesh_distance,
        symprec: float = defaults.symmetry_length_tolerance,
        angle_tolerance: float = defaults.symmetry_angle_tolerance,
        time_reversal: bool = True,
    ) -> None:
        """Initialize the generator.

        Args:
            initial_structure: Input crystal structure.
            task: Calculation task type determining defaults.
            kpt_density: K-point density in Å (reciprocal space).
            kpt_mode: K-point generation mode (None for task default).
            gamma_centered: Force gamma-centered mesh (VASP convention).
            only_even_num_kpts: Round k-point numbers up to even values.
            num_kpt_factor: K-point density multiplier (also sets NKRED).
            band_ref_dist: Reference distance for band path k-points.
            symprec: Symmetry tolerance for lengths (Å).
            angle_tolerance: Symmetry tolerance for angles (degrees).
            time_reversal: Whether system has time-reversal symmetry.
        """
        if kpt_mode:
            kpt_mode = KpointsMode.from_string(str(kpt_mode))

        self._initial_structure = initial_structure.copy()
        self._task = task
        self._kpt_density = kpt_density
        self._time_reversal = time_reversal
        self._symprec = symprec
        self._angle_tolerance = angle_tolerance

        # Set k-point factor with logging
        self._num_kpt_factor = num_kpt_factor or self._task.default_kpt_factor
        if self._num_kpt_factor != 1:
            logger.info(f"K-point factor is set to {self._num_kpt_factor}")

        # Initialize symmetry analyzer
        self._symmetrizer = StructureSymmetrizer(
            structure=self._initial_structure,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            band_mesh_distance=band_ref_dist,
            time_reversal=time_reversal,
        )

        # Apply task requirements (override user options if needed)
        self._gamma_centered = task.requisite_gamma_centered or gamma_centered
        self._adjust_only_even_num_kpts(only_even_num_kpts)
        self._adjust_kpt_mode(kpt_mode)

    def _adjust_only_even_num_kpts(self, user_option: bool) -> None:
        """Apply even k-point requirement, respecting task requisites."""
        task_requisite = self._task.requisite_only_even_num_kpts
        if task_requisite:
            logger.info("K-point numbers will be rounded to even values.")

        if isinstance(task_requisite, bool):
            self._only_even_num_kpts = task_requisite
        else:
            self._only_even_num_kpts = user_option

    def _adjust_kpt_mode(self, user_option: Optional[KpointsMode]) -> None:
        """Apply k-point mode, using task default if not specified."""
        self._kpt_mode = user_option or self._task.default_kpt_mode

    @property
    def _longest_axis_index(self) -> int:
        """Index of the longest lattice axis (for 2D defect handling)."""
        return int(np.argmax(self.structure.lattice.abc))

    def generate_input(self) -> None:
        """Generate both structure and k-points."""
        self._make_structure()
        self._make_kpoints()

    def _make_structure(self) -> None:
        """Prepare structure based on k-point mode."""
        if self._kpt_mode == KpointsMode.uniform:
            self._structure = self._initial_structure.copy()
        elif self._kpt_mode in (KpointsMode.primitive, KpointsMode.band):
            self._structure = self._symmetrizer.primitive
        else:
            raise NotImplementedError(f"Unknown k-point mode: {self._kpt_mode}")

        if self._structure != self._initial_structure:
            logger.info(
                "Input structure was converted to primitive cell. "
                "Use --uniform_kpt_mode to keep original structure."
            )

    def _make_kpoints(self) -> None:
        """Generate k-points based on settings."""
        self._set_num_kpt_list()
        self._set_kpt_shift()
        self._set_kpoints()

    def _set_num_kpt_list(self) -> None:
        """Calculate number of k-points along each axis."""
        # Use task-required values if specified
        if self._task.requisite_num_kpt_list:
            self._num_kpt_list = self._task.requisite_num_kpt_list
            return

        kpt_list: List[int] = []
        for reciprocal_length in self._reciprocal_lattice_lengths:
            raw_kpts = self._kpt_density * reciprocal_length

            if self._only_even_num_kpts:
                num_kpts = ceil(raw_kpts / 2) * 2
            else:
                num_kpts = ceil(raw_kpts)

            kpt_list.append(num_kpts * self._num_kpt_factor)

        # Special handling for 2D defects
        if self._task == Task.defect_2d:
            axis = "xyz"[self._longest_axis_index]
            logger.info(f"K-points along {axis} direction set to 1 for 2D defect.")
            kpt_list[self._longest_axis_index] = 1

        self._num_kpt_list = kpt_list

    @property
    def _reciprocal_lattice_lengths(self) -> Tuple[float, ...]:
        """Get reciprocal lattice vector lengths for k-point calculation.

        For oI and tI Bravais lattices, uses geometric mean to ensure
        equal k-points in all directions (required for symmetry).
        """
        reciprocal_abc = self._structure.lattice.reciprocal_lattice.abc
        self.bravais = BravaisLattice.from_sg_num(self._symmetrizer.sg_number)

        if self._kpt_mode.band_or_primitive and self.bravais.need_same_num_kpt:
            logger.warning(
                "For oI and tI Bravais lattices, equal k-points in all "
                "directions are required to preserve symmetry."
            )
            geometric_mean = pow(float(np.prod(reciprocal_abc)), 1 / 3)
            return (geometric_mean,) * 3

        return reciprocal_abc

    def _set_kpt_shift(self) -> None:
        """Determine k-point mesh shift."""
        if self._gamma_centered:
            self._kpt_shift = [0.0, 0.0, 0.0]
            return

        if self._kpt_mode is KpointsMode.uniform:
            # Shift along directions perpendicular to 90-degree angles
            kpt_shift = []
            angles = self._structure.lattice.angles
            for i in range(3):
                is_perpendicular = (
                    abs(angles[i - 2] - 90) < 1e-5
                    and abs(angles[i - 1] - 90) < 1e-5
                )
                kpt_shift.append(0.5 if is_perpendicular else 0.0)
        else:
            kpt_shift = self.bravais.kpt_centering

        # Force no shift along vacuum direction for 2D defects
        if (
            kpt_shift[self._longest_axis_index] != 0.0
            and self._task is Task.defect_2d
        ):
            kpt_shift[self._longest_axis_index] = 0.0

        # Only apply shift for even k-point counts
        even_kpts = [n % 2 == 0 for n in self._num_kpt_list]
        self._kpt_shift = [s * int(e) for s, e in zip(kpt_shift, even_kpts)]

    def _set_kpoints(self) -> None:
        """Generate final Kpoints object."""
        # Create fresh symmetrizer for possibly transformed structure
        symmetrizer = StructureSymmetrizer(
            self._structure,
            symprec=self._symprec,
            angle_tolerance=self._angle_tolerance,
            time_reversal=self._time_reversal,
        )
        irreducible_kpoints = symmetrizer.irreducible_kpoints(
            self._num_kpt_list, self._kpt_shift
        )

        self.comment = ""

        if self._kpt_mode is KpointsMode.band:
            self._add_weighted_kpts(irreducible_kpoints)
            self._add_band_path_kpts()
            self._num_kpts = len(self._kpoints.kpts)
        else:
            self._kpoints = Kpoints(
                comment=self.comment,
                kpts=(self._num_kpt_list,),
                kpts_shift=self._kpt_shift,
            )
            self._num_kpts = len(irreducible_kpoints)
            self.kpoints.comment += f"Num irrep kpoints: {self._num_kpts}"

    def _add_weighted_kpts(
        self, irreducible_kpoints: List[Tuple[List[float], int]]
    ) -> None:
        """Add weighted irreducible k-points for band calculations."""
        kpts = [k[0] for k in irreducible_kpoints]
        weights = [k[1] for k in irreducible_kpoints]
        labels = [None] * len(irreducible_kpoints)

        self._kpoints = Kpoints(
            comment=self.comment,
            style=Kpoints.supported_modes.Reciprocal,
            num_kpts=len(kpts),
            kpts=kpts,
            kpts_weights=weights,
            labels=labels,
        )

    def _add_band_path_kpts(self) -> None:
        """Append band structure path k-points from seekpath."""
        k_path = self._symmetrizer.seekpath_data["explicit_kpoints_rel"]
        k_labels = self._symmetrizer.seekpath_data["explicit_kpoints_labels"]
        formula = self._structure.composition.reduced_formula
        sg = self._symmetrizer.sg_number

        self.kpoints.comment += (
            f"k-path added by seekpath. Formula: {formula} SG: {sg} "
        )
        self._kpoints.num_kpts += len(k_path)
        self._kpoints.kpts += list(k_path)
        self._kpoints.labels += k_labels
        self._kpoints.kpts_weights += [0] * len(k_path)

    @property
    def structure(self) -> Structure:
        """The generated structure."""
        return self._structure

    @property
    def kpoints(self) -> Kpoints:
        """The generated Kpoints object."""
        return self._kpoints

    @property
    def num_kpts(self) -> int:
        """Number of irreducible k-points."""
        return self._num_kpts

    @property
    def num_kpt_factor(self) -> int:
        """K-point density multiplier."""
        return self._num_kpt_factor
