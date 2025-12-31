# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""VASP input file generation and management.

This module provides the VaspInputFiles class that coordinates the generation
of all VASP input files (INCAR, POSCAR, KPOINTS, POTCAR) from categorized
input options.

Key Components:
    - VaspInputFiles: Main class for generating complete VASP input sets
"""
from __future__ import annotations

from copy import copy
from pathlib import Path
from typing import Any, Dict, Optional

from pymatgen.io.vasp.sets import Kpoints, Poscar, Potcar

from vise import __version__
from vise.input_set.incar import ViseIncar
from vise.input_set.incar_settings_generator import IncarSettingsGenerator
from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.kpoints import ViseKpoints
from vise.input_set.potcar_generator import generate_potcar
from vise.input_set.structure_kpoints_generator import StructureKpointsGenerator
from vise.input_set.vise_log import ViseLog
from vise.util.logger import get_logger
from vise.util.structure_handler import create_symbol_list

logger = get_logger(__name__)


class VaspInputFiles:
    """Generator and container for complete VASP input file sets.

    This class coordinates the generation of all VASP input files by:
    1. Generating structure and k-points from StructureKpointsGenerator
    2. Generating POTCAR from potcar_generator
    3. Generating INCAR settings from IncarSettingsGenerator
    4. Applying any user-specified INCAR overrides

    Attributes:
        incar: The generated ViseIncar object.
        poscar: The generated Poscar object.
        kpoints: The generated Kpoints object.
        potcar: The generated Potcar object.
        vise_log: ViseLog object for provenance tracking.
        initial_structure: The original input structure.
        version: The vise version string.
    """

    def __init__(
        self,
        input_options: CategorizedInputOptions,
        overridden_incar_settings: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Initialize and generate VASP input files.

        Args:
            input_options: Categorized input options containing structure,
                task, XC functional, and other generation parameters.
            overridden_incar_settings: Optional dictionary of INCAR settings
                to override the automatically generated values.
        """
        self._version = __version__
        self._input_options = input_options
        self._initial_structure = input_options.initial_structure
        self._generate_structure_kpoints(input_options)
        self._generate_potcar_incar_settings(input_options)
        self._overridden_incar_settings = overridden_incar_settings
        self._incar_settings.update(overridden_incar_settings or {})

    def _generate_structure_kpoints(
        self,
        input_options: CategorizedInputOptions,
    ) -> None:
        """Generate structure and k-points from input options.

        Args:
            input_options: Categorized input options.
        """
        str_kpoints_generator = StructureKpointsGenerator(
            **input_options.structure_kpoints_options
        )
        str_kpoints_generator.generate_input()

        self._structure = str_kpoints_generator.structure
        self._kpoints = str_kpoints_generator.kpoints
        self._num_kpts = str_kpoints_generator.num_kpts
        self._num_kpt_factor = str_kpoints_generator.num_kpt_factor

    def _generate_potcar_incar_settings(
        self,
        input_options: CategorizedInputOptions,
    ) -> None:
        """Generate POTCAR and INCAR settings from input options.

        Args:
            input_options: Categorized input options.
        """
        symbol_list = create_symbol_list(self._structure)
        self._potcar = generate_potcar(
            symbol_list=symbol_list,
            **input_options.potcar_options,
        )

        incar_settings_generator = IncarSettingsGenerator(
            self._structure,
            symbol_list,
            self._num_kpts,
            self._num_kpt_factor,
            self._potcar,
            **input_options.incar_settings_options,
        )
        self._incar_settings = incar_settings_generator.incar_settings

    @property
    def vise_log(self) -> ViseLog:
        """Return ViseLog with user INCAR settings added.

        Returns:
            ViseLog object for provenance tracking.
        """
        result = copy(self._input_options.vise_log)
        result.user_incar_settings = self._overridden_incar_settings
        return result

    @property
    def initial_structure(self):
        """Return the original input structure.

        Returns:
            The initial pymatgen Structure object.
        """
        return self._initial_structure

    @property
    def version(self) -> str:
        """Return the vise version string.

        Returns:
            Version string.
        """
        return self._version

    @property
    def incar(self) -> ViseIncar:
        """Return the generated ViseIncar object.

        Returns:
            ViseIncar object with all settings applied.
        """
        return ViseIncar.from_dict(self._incar_settings)

    @property
    def poscar(self) -> Poscar:
        """Return the generated Poscar object.

        Returns:
            pymatgen Poscar object.
        """
        return Poscar(self._structure)

    @property
    def kpoints(self) -> Kpoints:
        """Return the generated Kpoints object.

        Returns:
            pymatgen Kpoints object.
        """
        return self._kpoints

    @property
    def potcar(self) -> Potcar:
        """Return the generated Potcar object.

        Returns:
            pymatgen Potcar object.
        """
        return self._potcar

    @property
    def input_files(self) -> Dict[str, Any]:
        """Return all input file objects as a dictionary.

        Returns:
            Dictionary mapping filenames to their objects.
        """
        return {
            "INCAR": self.incar,
            "KPOINTS": self.kpoints,
            "POSCAR": self.poscar,
            "POTCAR": self.potcar,
        }

    def create_input_files(
        self,
        dirname: Path,
        poscar_significant_figures: int = 10,
    ) -> None:
        """Write all VASP input files to a directory.

        Args:
            dirname: Directory path to write files to (created if needed).
            poscar_significant_figures: Number of significant figures for
                atomic coordinates in POSCAR (default: 10).
        """
        dirname.mkdir(exist_ok=True)

        for filename, obj in self.input_files.items():
            if isinstance(obj, Poscar):
                # Write POSCAR with higher precision than pymatgen default
                obj.write_file(
                    dirname / filename,
                    significant_figures=poscar_significant_figures,
                )
            elif isinstance(obj, Kpoints):
                # Use ViseKpoints for improved formatting
                kpoints = ViseKpoints.from_dict(obj.as_dict())
                kpoints.write_file(dirname / filename)
            else:
                obj.write_file(dirname / filename)
