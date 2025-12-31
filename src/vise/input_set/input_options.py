# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Input options management and categorization for VASP input generation.

This module provides utilities for managing and categorizing input options
that are distributed across different VASP input generation components
(structure/k-points, POTCAR, INCAR settings).

Key Components:
    - CategorizedInputOptions: Categorizes and distributes input options
    - ViseInputOptionsError: Exception for invalid input options
"""
from __future__ import annotations

from copy import copy, deepcopy
from inspect import getfullargspec
from typing import Any, Dict, List, Set

from pymatgen.core.structure import Structure

from vise import __version__
from vise.analyzer.band_edge_properties import is_band_gap
from vise.defaults import defaults
from vise.input_set.incar_settings_generator import IncarSettingsGenerator
from vise.input_set.potcar_generator import generate_potcar
from vise.input_set.structure_kpoints_generator import StructureKpointsGenerator
from vise.input_set.task import Task
from vise.input_set.vise_log import ViseLog
from vise.input_set.xc import Xc
from vise.util.logger import get_logger

logger = get_logger(__name__)


def function_args(func) -> List[str]:
    """Extract argument names from a function.

    Args:
        func: A function object.

    Returns:
        List of argument names.
    """
    return getfullargspec(func).args


def constructor_args(target_class) -> List[str]:
    """Extract constructor argument names from a class, excluding 'self'.

    Args:
        target_class: A class object.

    Returns:
        List of constructor argument names (excluding 'self').
    """
    result = function_args(target_class.__init__)
    result.remove("self")
    return result


# Extract argument names from each input generation component
potcar_args: List[str] = function_args(generate_potcar)
incar_settings_args: List[str] = constructor_args(IncarSettingsGenerator)
structure_kpoints_args: List[str] = constructor_args(StructureKpointsGenerator)

#: Set of all valid input option names
assignable_option_set: Set[str] = (
    set(structure_kpoints_args) | set(potcar_args) | set(incar_settings_args)
)


class CategorizedInputOptions:
    """Categorize and distribute input options to appropriate generators.

    This class takes a flat dictionary of input options and categorizes them
    into groups that correspond to different VASP input generation components:
    - Structure and k-points generation options
    - POTCAR generation options
    - INCAR settings generation options

    It also handles automatic k-point density selection based on band gap
    information and task type.

    Attributes:
        initial_structure: The original input structure (copied).
        task: The calculation task type.
        xc: The exchange-correlation functional.
    """

    def __init__(
        self,
        structure: Structure,
        task: Task,
        xc: Xc,
        **input_options: Any,
    ) -> None:
        """Initialize with structure, task, XC functional, and options.

        Args:
            structure: pymatgen Structure object.
            task: Calculation task type.
            xc: Exchange-correlation functional.
            **input_options: Additional input options to be categorized.

        Raises:
            ViseInputOptionsError: If unknown options are provided.
        """
        self._input_options = deepcopy(input_options)
        self.initial_structure = structure.copy()
        self.task = task
        self.xc = xc
        self._set_kpt_density()
        self._raise_error_when_unknown_options_exist()

    def _raise_error_when_unknown_options_exist(self) -> None:
        """Validate that all provided options are recognized.

        Raises:
            ViseInputOptionsError: If any unknown options are found.
        """
        unknown_args_set = self.input_option_set - assignable_option_set
        if unknown_args_set:
            raise ViseInputOptionsError(f"Options {unknown_args_set} are invalid")

    def _set_kpt_density(self) -> None:
        """Set k-point density based on band gap and task type.

        Automatically determines appropriate k-point density:
        - Lower density for defect calculations
        - Higher density for insulators (has band gap)
        - Default density for metals
        """
        band_gap = self._input_options.get("band_gap", None)
        vbm_cbm = self._input_options.get("vbm_cbm", None)
        kpt_density = self._input_options.get("kpt_density", None)

        msg = ""
        if kpt_density is None:
            if self.task in (Task.defect, Task.defect_2d):
                kpt_density = defaults.defect_kpoint_density
                msg = " because task is defect"
            elif is_band_gap(band_gap, vbm_cbm):
                kpt_density = defaults.insulator_kpoint_density
                msg = " because there is a band gap"
            else:
                kpt_density = defaults.kpoint_density

            self._input_options["kpt_density"] = kpt_density

        logger.info(f"K-point density is set to {kpt_density}{msg}.")

    @property
    def all_options(self) -> Dict[str, Any]:
        """Return all options including structure, task, and XC functional.

        Returns:
            Dictionary containing all options and core parameters.
        """
        result = copy(self._input_options)
        result.update({
            "task": self.task,
            "xc": self.xc,
            "initial_structure": self.initial_structure,
        })
        return result

    @property
    def input_option_set(self) -> Set[str]:
        """Return the set of provided input option names.

        Returns:
            Set of option names that were provided.
        """
        return set(self._input_options.keys())

    def pick_target_options(self, target_args: List[str]) -> Dict[str, Any]:
        """Extract options that match the target argument names.

        Args:
            target_args: List of argument names to extract.

        Returns:
            Dictionary of matching options.
        """
        result = {}
        for target_arg in target_args:
            if target_arg in self.all_options:
                result[target_arg] = self.all_options[target_arg]
        return result

    @property
    def structure_kpoints_options(self) -> Dict[str, Any]:
        """Return options for StructureKpointsGenerator.

        Returns:
            Dictionary of structure/k-points generation options.
        """
        return self.pick_target_options(target_args=structure_kpoints_args)

    @property
    def potcar_options(self) -> Dict[str, Any]:
        """Return options for POTCAR generation.

        Returns:
            Dictionary of POTCAR generation options.
        """
        return self.pick_target_options(target_args=potcar_args)

    @property
    def incar_settings_options(self) -> Dict[str, Any]:
        """Return options for IncarSettingsGenerator.

        Returns:
            Dictionary of INCAR settings generation options.
        """
        return self.pick_target_options(target_args=incar_settings_args)

    @property
    def vise_log(self) -> ViseLog:
        """Create a ViseLog object with current options.

        Returns:
            ViseLog object for provenance tracking.
        """
        return ViseLog(
            version=__version__,
            task=self.task,
            xc=self.xc,
            input_options=self._input_options,
        )


class ViseInputOptionsError(KeyError):
    """Exception raised for invalid input options."""

    pass
