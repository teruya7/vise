# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Prior calculation information for VASP.

This module provides PriorInfo for storing and using information from
previous calculations to configure new VASP runs.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml
from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen.core import Structure
from pymatgen.io.vasp import Outcar, Potcar, Vasprun

from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.defaults import defaults


@dataclass
class PriorInfo(MSONable):
    """Previous calculation information for guiding new calculations.

    Stores results from prior calculations that can be used to
    improve settings for subsequent VASP runs.

    Attributes:
        structure: Crystal structure from previous calculation.
        energy_per_atom: Total energy per atom (eV).
        band_gap: Band gap in eV.
        vbm_cbm: VBM and CBM energies as [vbm, cbm].
        total_magnetization: Total magnetic moment in μB.
        data_source: Source of the data (e.g., "vasp", "mp").
        is_cluster: Whether structure is an isolated molecule/cluster.
        charge: Net charge on the system.
        icsd_ids: ICSD database IDs if applicable.
        incar: INCAR settings to inherit.
    """

    structure: Optional[Structure] = None
    energy_per_atom: Optional[float] = None
    band_gap: Optional[float] = None
    vbm_cbm: List[float] = field(default_factory=list)
    total_magnetization: Optional[float] = None
    data_source: Optional[str] = None
    is_cluster: Optional[bool] = None
    charge: Optional[int] = None
    icsd_ids: Optional[List[int]] = None
    incar: Dict[str, Any] = field(default_factory=dict)

    def dump_yaml(self, filename: str = "prior_info.yaml") -> None:
        """Save to YAML file."""
        with open(filename, "w") as f:
            f.write(yaml.dump(self.as_dict()))

    @classmethod
    def load_yaml(cls, filename: str = "prior_info.yaml") -> "PriorInfo":
        """Load from YAML file."""
        with open(filename, "r") as f:
            data = yaml.load(f, Loader=yaml.SafeLoader)
        return cls.from_dict(data)

    def dump_json(self, filename: str = "prior_info.json") -> None:
        """Save to JSON file."""
        with open(filename, "w") as f:
            json.dump(self.as_dict(), f, indent=2)

    @classmethod
    def load_json(cls, filename: str = "prior_info.json") -> "PriorInfo":
        """Load from JSON file."""
        return loadfn(filename)

    @property
    def is_magnetic(self) -> Optional[bool]:
        """Check if system has non-zero magnetization."""
        if self.total_magnetization is None:
            return None
        return abs(self.total_magnetization) > defaults.integer_criterion

    @property
    def has_band_gap(self) -> bool:
        """Check if system has a band gap above threshold."""
        return self.band_gap > defaults.band_gap_criterion

    @property
    def is_metal(self) -> bool:
        """Check if system is metallic (no band gap)."""
        return not self.has_band_gap

    @property
    def input_options_kwargs(self) -> Dict[str, Any]:
        """Get input options derived from prior info.

        Returns:
            Dictionary of options suitable for input generation.
        """
        result: Dict[str, Any] = {}

        if self.vbm_cbm:
            result["vbm_cbm"] = self.vbm_cbm
        if isinstance(self.is_magnetic, bool):
            result["time_reversal"] = not self.is_magnetic
        if self.band_gap:
            result["band_gap"] = self.band_gap
        if self.charge:
            result["charge"] = self.charge

        return result


def prior_info_from_calc_dir(
    prev_dir_path: Path,
    vasprun: str = "vasprun.xml",
    outcar: str = "OUTCAR",
    potcar: str = "POTCAR",
) -> PriorInfo:
    """Extract PriorInfo from a completed VASP calculation directory.

    Args:
        prev_dir_path: Path to the calculation directory.
        vasprun: Name of the vasprun.xml file.
        outcar: Name of the OUTCAR file.
        potcar: Name of the POTCAR file.

    Returns:
        PriorInfo populated with calculation results.
    """
    vasprun_obj = Vasprun(str(prev_dir_path / vasprun))
    outcar_obj = Outcar(str(prev_dir_path / outcar))
    potcar_obj = Potcar.from_file(str(prev_dir_path / potcar))

    charge = get_net_charge_from_vasp(
        vasprun_obj.final_structure,
        vasprun_obj.parameters["NELECT"],
        potcar_obj,
    )

    structure = vasprun_obj.final_structure.copy()
    energy_per_atom = outcar_obj.final_energy / len(structure)
    band_edge_property = VaspBandEdgeProperties(vasprun_obj, outcar_obj)

    return PriorInfo(
        structure=structure,
        charge=charge,
        energy_per_atom=energy_per_atom,
        band_gap=band_edge_property.band_gap,
        vbm_cbm=band_edge_property.vbm_cbm,
        total_magnetization=outcar_obj.total_mag,
    )


def get_net_charge_from_vasp(
    structure: Structure,
    nelect: int,
    potcar: Potcar,
) -> int:
    """Calculate net charge from VASP output.

    Compares the number of electrons in NELECT with the nominal
    valence electron count from POTCARs.

    Args:
        structure: Crystal structure.
        nelect: Number of electrons from INCAR/OUTCAR.
        potcar: POTCAR object with pseudopotential info.

    Returns:
        Net charge (positive = electron deficit).

    Raises:
        ValueError: If POTCAR elements don't match structure.
    """
    nuclei_charge = 0

    for elem, pot in zip(structure.composition, potcar):
        if pot.element != str(elem):
            raise ValueError(
                f"Element mismatch: POTCAR has {pot.element}, "
                f"structure has {elem}"
            )
        nuclei_charge += pot.nelectrons * structure.composition[elem]

    # Charge is difference between nuclear and electronic charge
    return int(nuclei_charge - nelect)
