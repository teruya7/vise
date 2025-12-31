# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""VASP calculation task type definitions.

This module defines the Task enum representing different types of VASP
calculations, each with specific settings for k-points, relaxation behavior,
and other parameters.
"""

from typing import List, Optional

from monty.json import MSONable

from vise.input_set.kpoints_mode import KpointsMode
from vise.util.enum import ExtendedEnum


class Task(MSONable, ExtendedEnum):
    """VASP calculation task types.

    Each task type determines appropriate INCAR settings, k-point
    configuration, and relaxation behavior.

    Attributes:
        structure_opt: Standard structure optimization
        structure_opt_rough: Quick structure optimization (looser convergence)
        structure_opt_tight: Tight structure optimization (strict convergence)
        cluster_opt: Isolated cluster/molecule optimization
        phonon_force: Force calculation for phonons
        defect: Point defect calculation
        defect_2d: Point defect in 2D material
        band: Band structure calculation
        dos: Density of states calculation
        dielectric_dfpt: Dielectric constant via DFPT
        dielectric_finite_field: Dielectric via finite field
        dielectric_function: Frequency-dependent dielectric function
    """

    structure_opt = "structure_opt"
    structure_opt_rough = "structure_opt_rough"
    structure_opt_tight = "structure_opt_tight"
    cluster_opt = "cluster_opt"
    phonon_force = "phonon_force"
    defect = "defect"
    defect_2d = "defect_2d"
    band = "band"
    dos = "dos"
    dielectric_dfpt = "dielectric_dfpt"
    dielectric_finite_field = "dielectric_finite_field"
    dielectric_function = "dielectric_function"

    @property
    def is_lattice_relaxed(self) -> bool:
        """Check if task involves lattice relaxation."""
        return "structure_opt" in self.name

    @property
    def is_atom_relaxed_lattice_fixed(self) -> bool:
        """Check if task relaxes atoms but keeps lattice fixed."""
        return self in (self.cluster_opt, self.defect, self.defect_2d)

    @property
    def is_atom_relaxed(self) -> bool:
        """Check if task involves any atomic relaxation."""
        return self.is_lattice_relaxed or self.is_atom_relaxed_lattice_fixed

    @property
    def is_dielectric(self) -> bool:
        """Check if task is a dielectric calculation."""
        return "dielectric" in self.name

    @property
    def is_tight_calc(self) -> bool:
        """Check if task requires tight convergence."""
        return self in (self.structure_opt_tight, self.phonon_force)

    @property
    def is_plot_task(self) -> bool:
        """Check if task produces plottable results."""
        return self in (self.band, self.dos, self.dielectric_function)

    @property
    def is_spectrum_task(self) -> bool:
        """Check if task produces spectrum data (DOS or dielectric)."""
        return self in (self.dos, self.dielectric_function)

    @property
    def default_kpt_factor(self) -> Optional[int]:
        """Get default k-point density multiplier for this task.

        Returns:
            K-point factor: 3 for dielectric_function, 2 for DOS/dielectric,
            1 for band/defect/optimization tasks.

        Raises:
            NotImplementedError: If factor not defined for this task.
        """
        if self is self.dielectric_function:
            return 3
        elif self in (self.dos, self.dielectric_dfpt, self.dielectric_finite_field):
            return 2
        elif self._condition_kpt_factor_is_one:
            return 1
        else:
            raise NotImplementedError(f"K-point factor not defined for {self}")

    @property
    def _condition_kpt_factor_is_one(self) -> bool:
        """Check if task uses factor 1 for k-points."""
        return (
            self in (self.band, self.cluster_opt, self.phonon_force,
                     self.defect, self.defect_2d)
            or self.is_lattice_relaxed
        )

    @property
    def default_kpt_mode(self) -> KpointsMode:
        """Get default k-point generation mode for this task.

        Returns:
            KpointsMode appropriate for this calculation type.

        Raises:
            NotImplementedError: If mode not defined for this task.
        """
        if self == self.band:
            return KpointsMode.band
        elif self in (self.defect, self.defect_2d, self.cluster_opt,
                      self.phonon_force):
            return KpointsMode.uniform
        elif self._condition_kpoints_mode_is_primitive:
            return KpointsMode.primitive
        else:
            raise NotImplementedError(f"K-point mode not defined for {self}")

    @property
    def _condition_kpoints_mode_is_primitive(self) -> bool:
        """Check if task uses primitive k-point mode."""
        return self is self.dos or self.is_dielectric or self.is_lattice_relaxed

    @property
    def requisite_num_kpt_list(self) -> Optional[List[int]]:
        """Get required k-point mesh, if task has a fixed requirement."""
        if self == self.cluster_opt:
            return [1, 1, 1]
        return None

    @property
    def requisite_only_even_num_kpts(self) -> Optional[bool]:
        """Check if task requires only even k-point numbers."""
        if self == self.cluster_opt:
            return False
        elif self == self.band:
            return True
        return None

    @property
    def requisite_gamma_centered(self) -> Optional[bool]:
        """Check if task requires gamma-centered k-mesh.

        Gamma-center is required for GW calculations and strongly
        recommended for DOS and dielectric function to sample band edges.
        """
        if self in (self.cluster_opt, self.phonon_force) or self.is_spectrum_task:
            return True
        return None

    @property
    def need_spin(self) -> bool:
        """Check if task typically requires spin polarization."""
        return self in (self.defect, self.defect_2d)

    @property
    def fine_to_inherit_structure_from_prev(self) -> bool:
        """Check if structure can be inherited from previous calculation.

        For phonon calculations, POSCAR must not be inherited to ensure
        correct displacement setup.
        """
        return self is not self.phonon_force


# Backward compatibility
condition_kpt_factor_is_one = Task._condition_kpt_factor_is_one
condition_kpoints_mode_is_primitive = Task._condition_kpoints_mode_is_primitive
