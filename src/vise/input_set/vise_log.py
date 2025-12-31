# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
"""VISE calculation log file handling.

This module provides the ViseLog dataclass for logging and persisting
VASP calculation settings, enabling reproducibility and tracking of
input parameters.
"""

from copy import copy
from dataclasses import dataclass
from typing import Any, Dict, Optional

import yaml
from monty.json import MSONable

from vise.input_set.datasets.potcar_set import PotcarSet
from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.mix_in import ToYamlFileMixIn


@dataclass
class ViseLog(MSONable, ToYamlFileMixIn):
    """Log of VISE calculation settings.

    Stores the parameters used to generate VASP input files for
    reproducibility and provenance tracking.

    Attributes:
        version: VISE version used for generation.
        task: Calculation task type.
        xc: Exchange-correlation functional.
        input_options: Dictionary of input generation options.
        user_incar_settings: User-specified INCAR overrides.
        ldauu: DFT+U parameters (U values by element).
        ldaul: DFT+U parameters (L values by element).
    """

    version: str
    task: Task
    xc: Xc
    input_options: Dict[str, Any]
    user_incar_settings: Optional[Dict[str, Any]] = None
    ldauu: Optional[Dict[str, float]] = None
    ldaul: Optional[Dict[str, float]] = None

    def to_yaml(self) -> str:
        """Serialize to YAML string with preserved key order."""
        return yaml.safe_dump(self.as_dict(), sort_keys=False)

    def as_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        result: Dict[str, Any] = {
            "version": self.version,
            "task": str(self.task),
            "xc": str(self.xc),
            "input_options": self.input_options,
        }

        if self.user_incar_settings is not None:
            result["user_incar_settings"] = self.user_incar_settings
        if self.ldauu is not None:
            result["ldauu"] = self.ldauu
        if self.ldaul is not None:
            result["ldaul"] = self.ldaul

        # Convert enum values to strings in input_options
        if "kpt_mode" in self.input_options:
            self.input_options["kpt_mode"] = str(self.input_options["kpt_mode"])
        if "potcar_set" in self.input_options:
            self.input_options["potcar_set"] = str(self.input_options["potcar_set"])

        return result

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ViseLog":
        """Create ViseLog from dictionary.

        Args:
            d: Dictionary with serialized ViseLog data.

        Returns:
            Reconstructed ViseLog instance.
        """
        data = copy(d)
        data["task"] = Task.from_string(d["task"])
        data["xc"] = Xc.from_string(d["xc"])

        # Convert string enums back to enum types
        input_options = data["input_options"]
        if "kpt_mode" in input_options:
            input_options["kpt_mode"] = KpointsMode.from_string(
                input_options["kpt_mode"]
            )
        if "potcar_set" in input_options:
            input_options["potcar_set"] = PotcarSet.from_string(
                input_options["potcar_set"]
            )

        return super().from_dict(data)
