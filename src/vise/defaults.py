# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Default configuration and settings management for VISE.

This module provides the singleton `Defaults` class that manages all default
values used throughout the VISE package. Defaults can be overridden via
user settings files (vise.yaml or .vise.yaml) found in the directory hierarchy.

Usage:
    from vise.defaults import defaults
    print(defaults.kpoint_density)  # Access default k-point density

Key Components:
    - DefaultsBase: Abstract base class for defaults management
    - Defaults: Singleton class containing all package defaults
    - defaults: Global singleton instance
"""
from __future__ import annotations

import warnings
from abc import ABC
from pathlib import Path
from typing import Any, Dict, List

from monty.design_patterns import singleton

from vise.input_set.datasets.potcar_set import PotcarSet
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.user_settings import UserSettings
from vise.util.logger import get_logger

logger = get_logger(__name__)


class DefaultsBase(ABC):
    """Abstract base class for default settings management.

    This class provides the mechanism for loading and applying user settings
    from YAML configuration files.

    Attributes:
        yaml_files: List of YAML files found in the directory hierarchy.
        user_settings: Dictionary of user-defined settings.
    """

    yaml_files: List[Path]
    user_settings: Dict[str, Any]

    def set_user_settings(self, yaml_filename: str) -> None:
        """Load user settings from YAML files in the directory hierarchy.

        Searches for YAML files (both with and without leading dot) from the
        current directory up to the root, and applies any matching settings
        to override default values.

        Args:
            yaml_filename: Name of the YAML settings file to search for
                (e.g., 'vise.yaml').
        """
        user_settings = UserSettings(yaml_filename=yaml_filename)
        self.yaml_files = user_settings.yaml_files_from_root_dir
        self.user_settings = user_settings.user_settings

        for key, value in self.user_settings.items():
            if hasattr(self, key):
                self.__setattr__("_" + key, value)
            else:
                logger.warning(f"{key} key in {yaml_filename} is not valid.")


@singleton
class Defaults(DefaultsBase):
    """Singleton class containing all VISE package default values.

    This class holds default values for symmetry tolerances, k-point densities,
    VASP file paths, POTCAR settings, and other configuration options. Values
    can be overridden via vise.yaml files in the directory hierarchy.

    The singleton pattern ensures consistent defaults across all VISE modules.

    Attributes:
        symmetry_length_tolerance: Tolerance for symmetry finding in Angstrom.
        symmetry_angle_tolerance: Tolerance for symmetry finding in degrees.
        dos_step_size: Energy step size for DOS calculations in eV.
        kpoint_density: Default k-point density for metals.
        insulator_kpoint_density: K-point density for insulators.
        defect_kpoint_density: K-point density for defect calculations.
        band_mesh_distance: Distance between k-points for band structure.
        str_opt_encut_factor: ENCUT multiplier for structure optimization.
        band_gap_criterion: Threshold for band gap detection in eV.
        integer_criterion: Tolerance for integer detection.
        default_num_cores: Default number of CPU cores for KPAR calculation.
        unused_core_ratio_threshold: Threshold for unused cores in KPAR.
        task: Default calculation task type.
        xc: Default exchange-correlation functional.
        options: User-defined additional options.
        user_incar_settings: User-defined INCAR overrides.
        user_incar_tags: User-defined INCAR tag categories.
        u_parameter_set_yaml_file: Path to Hubbard U parameter file.
        outcar: Path to OUTCAR file.
        contcar: Path to CONTCAR file.
        vasprun: Path to vasprun.xml file.
        procar: Path to PROCAR file.
        overridden_potcar: List of POTCAR overrides.
        potcar_set: POTCAR selection preset.
        lda_potcar: LDA POTCAR directory name.
        gga_potcar: GGA POTCAR directory name.
    """

    def __init__(self) -> None:
        """Initialize default values and load user settings."""
        # Symmetry detection tolerances
        self._symmetry_length_tolerance: float = 0.01
        self._symmetry_angle_tolerance: float = 5.0

        # K-point and energy grid settings
        self._dos_step_size: float = 0.01
        self._kpoint_density: float = 5.0
        self._insulator_kpoint_density: float = 2.5
        self._defect_kpoint_density: float = 1.8
        self._band_mesh_distance: float = 0.025

        # Calculation parameters
        self._str_opt_encut_factor: float = 1.3
        self._band_gap_criterion: float = 0.2  # in eV
        self._integer_criterion: float = 0.1
        self._default_num_cores: int = 4
        self._unused_core_ratio_threshold: float = 0.25

        # Default task and functional
        self._task: str = str(Task.structure_opt)
        self._xc: str = str(Xc.pbe)

        # User customization options
        self._options: Dict[str, Any] = {}
        self._user_incar_settings: Dict[str, Any] = {}
        self._user_incar_tags: Dict[str, List[str]] = {}

        # File paths
        self._u_parameter_set_yaml_file: Path = (
            Path(__file__).parent / "input_set" / "datasets" / "u_parameter_set.yaml"
        )
        self._outcar: str = "OUTCAR"
        self._contcar: str = "CONTCAR"
        self._vasprun: str = "vasprun.xml"
        self._procar: str = "PROCAR"

        # POTCAR settings
        self._overridden_potcar: List[str] = []
        self._potcar_set: str = str(PotcarSet.normal)
        self._lda_potcar: str = "LDA"
        self._gga_potcar: str = "PBE_54"

        # Warning suppression
        self._suppress_user_warning: bool = True

        # Load user settings from vise.yaml files
        self.set_user_settings(yaml_filename="vise.yaml")

        if self._suppress_user_warning:
            warnings.simplefilter("ignore", UserWarning)

    @property
    def symmetry_length_tolerance(self) -> float:
        """Tolerance for symmetry finding in Angstroms."""
        return self._symmetry_length_tolerance

    @property
    def symmetry_angle_tolerance(self) -> float:
        """Tolerance for symmetry finding in degrees."""
        return self._symmetry_angle_tolerance

    @property
    def dos_step_size(self) -> float:
        """Energy step size for DOS calculations in eV."""
        return self._dos_step_size

    @property
    def kpoint_density(self) -> float:
        """Default k-point density for metallic systems."""
        return self._kpoint_density

    @property
    def insulator_kpoint_density(self) -> float:
        """K-point density for insulating systems."""
        return self._insulator_kpoint_density

    @property
    def defect_kpoint_density(self) -> float:
        """K-point density for defect calculations."""
        return self._defect_kpoint_density

    @property
    def band_mesh_distance(self) -> float:
        """Distance between k-points for band structure calculations."""
        return self._band_mesh_distance

    @property
    def str_opt_encut_factor(self) -> float:
        """ENCUT multiplier for structure optimization."""
        return self._str_opt_encut_factor

    @property
    def band_gap_criterion(self) -> float:
        """Threshold for band gap detection in eV."""
        return self._band_gap_criterion

    @property
    def integer_criterion(self) -> float:
        """Tolerance for detecting integer values."""
        return self._integer_criterion

    @property
    def default_num_cores(self) -> int:
        """Default number of CPU cores for KPAR calculation."""
        return self._default_num_cores

    @property
    def unused_core_ratio_threshold(self) -> float:
        """Maximum acceptable ratio of unused cores for KPAR."""
        return self._unused_core_ratio_threshold

    @property
    def task(self) -> Task:
        """Default calculation task type."""
        return Task(self._task)

    @property
    def xc(self) -> Xc:
        """Default exchange-correlation functional."""
        return Xc(self._xc)

    @property
    def options(self) -> Dict[str, Any]:
        """User-defined additional options."""
        return self._options

    @property
    def user_incar_settings(self) -> Dict[str, Any]:
        """User-defined INCAR setting overrides."""
        return self._user_incar_settings

    @property
    def user_incar_tags(self) -> Dict[str, List[str]]:
        """User-defined INCAR tag categories."""
        return self._user_incar_tags

    @property
    def u_parameter_set_yaml_file(self) -> Path:
        """Path to the Hubbard U parameter configuration file."""
        return self._u_parameter_set_yaml_file

    @property
    def outcar(self) -> Path:
        """Path to the OUTCAR file."""
        return Path(self._outcar)

    @property
    def contcar(self) -> Path:
        """Path to the CONTCAR file."""
        return Path(self._contcar)

    @property
    def vasprun(self) -> Path:
        """Path to the vasprun.xml file."""
        return Path(self._vasprun)

    @property
    def procar(self) -> Path:
        """Path to the PROCAR file."""
        return Path(self._procar)

    @property
    def overridden_potcar(self) -> List[str]:
        """List of element-specific POTCAR overrides."""
        return self._overridden_potcar

    @property
    def potcar_set(self) -> PotcarSet:
        """POTCAR selection preset."""
        return PotcarSet(self._potcar_set)

    @property
    def lda_potcar(self) -> str:
        """LDA POTCAR directory name."""
        return self._lda_potcar

    @property
    def gga_potcar(self) -> str:
        """GGA POTCAR directory name."""
        return self._gga_potcar


#: Global singleton instance of Defaults
defaults: Defaults = Defaults()
