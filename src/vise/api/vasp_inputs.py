# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""VASP input generation API.

This module provides functions for generating VASP input files
without CLI dependencies.

Example:
    >>> from pymatgen.core import Structure
    >>> from vise.api.vasp_inputs import create_vasp_set
    >>>
    >>> struct = Structure.from_file("POSCAR")
    >>> inputs = create_vasp_set(struct, task="structure_opt", xc="pbe")
    >>> inputs.write("./calculation")
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any, Union

from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints, Poscar, Potcar

from vise.input_set.incar import ViseIncar
from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.kpoints import ViseKpoints
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.input_set.vasp_input_files import VaspInputFiles


@dataclass
class VaspInputSet:
    """Container for VASP input files.
    
    Attributes:
        incar: INCAR settings as ViseIncar object
        kpoints: KPOINTS as Kpoints object
        poscar: POSCAR as Poscar object
        potcar: POTCAR as Potcar object
        structure: Final structure used for input generation
        initial_structure: Original input structure
    """
    incar: ViseIncar
    kpoints: Kpoints
    poscar: Poscar
    potcar: Potcar
    structure: Structure
    initial_structure: Structure
    
    def write(
        self, 
        directory: Union[str, Path],
        poscar_significant_figures: int = 10
    ) -> None:
        """
        Write all VASP input files to a directory.
        
        Args:
            directory: Target directory path (created if doesn't exist)
            poscar_significant_figures: Number of significant figures for 
                POSCAR coordinates
        
        Example:
            >>> inputs = create_vasp_set(structure, task="structure_opt")
            >>> inputs.write("./my_calculation")
            # Creates: ./my_calculation/INCAR, KPOINTS, POSCAR, POTCAR
        """
        directory = Path(directory)
        directory.mkdir(exist_ok=True, parents=True)
        
        self.incar.write_file(directory / "INCAR")
        
        kpoints = ViseKpoints.from_dict(self.kpoints.as_dict())
        kpoints.write_file(directory / "KPOINTS")
        
        self.poscar.write_file(
            directory / "POSCAR", 
            significant_figures=poscar_significant_figures
        )
        self.potcar.write_file(directory / "POTCAR")
    
    def as_dict(self) -> Dict[str, Any]:
        """
        Get input files as dictionary for inspection.
        
        Returns:
            Dictionary with INCAR settings, KPOINTS info, etc.
        """
        return {
            "incar": dict(self.incar),
            "kpoints": self.kpoints.as_dict(),
            "structure": self.structure.as_dict(),
            "num_atoms": len(self.structure),
        }


def create_vasp_set(
    structure: Structure,
    task: Union[str, Task] = "structure_opt",
    xc: Union[str, Xc] = "pbe",
    kpt_density: Optional[float] = None,
    overridden_potcar: Optional[Dict[str, str]] = None,
    user_incar_settings: Optional[Dict[str, Any]] = None,
    **options
) -> VaspInputSet:
    """
    Create VASP input files for a given structure and task.
    
    This is the main entry point for programmatic VASP input generation.
    
    Args:
        structure: pymatgen Structure object
        task: Task type. Options:
            - "structure_opt" (default): Structure optimization
            - "band": Band structure calculation
            - "dos": Density of states
            - "dielectric_dfpt": Dielectric function (DFPT)
            - "dielectric_finite_field": Dielectric (finite field)
            - "phonon_force": Phonon force calculation
            - "defect": Defect calculation
            - "cluster_opt": Cluster optimization
        xc: Exchange-correlation functional. Options:
            - "pbe" (default): PBE-GGA
            - "hse": HSE06 hybrid
            - "scan": SCAN meta-GGA
            - "pbesol": PBEsol
            - "lda": LDA
        kpt_density: K-point density in Å. If None, automatically
            determined based on task and band gap.
        overridden_potcar: Override POTCAR for specific elements.
            Example: {"Mg": "Mg_pv", "O": "O_h"}
        user_incar_settings: Additional INCAR settings to override defaults.
            Example: {"EDIFF": 1e-6, "NELM": 200}
        **options: Additional options passed to CategorizedInputOptions.
            Common options include:
            - band_gap: Band gap in eV (affects k-point density)
            - charge: System charge
            - set_spin_orbit: Enable spin-orbit coupling
    
    Returns:
        VaspInputSet containing all input files
    
    Example:
        >>> from pymatgen.core import Structure
        >>> from vise.api.vasp_inputs import create_vasp_set
        >>>
        >>> # Basic structure optimization
        >>> struct = Structure.from_file("POSCAR")
        >>> inputs = create_vasp_set(struct)
        >>> inputs.write("./relax")
        >>>
        >>> # HSE band structure with custom settings
        >>> inputs = create_vasp_set(
        ...     struct,
        ...     task="band",
        ...     xc="hse",
        ...     user_incar_settings={"NBANDS": 100}
        ... )
        >>> inputs.write("./band")
        >>>
        >>> # Inspect INCAR settings before writing
        >>> print(inputs.incar)
    
    Raises:
        ValueError: If invalid task or xc is provided
        ViseInputOptionsError: If unknown options are provided
    """
    # Convert string to enum if needed
    if isinstance(task, str):
        task = Task(task)
    if isinstance(xc, str):
        xc = Xc(xc)
    
    # Set k-point density if provided
    if kpt_density is not None:
        options["kpt_density"] = kpt_density
    
    # Create input options
    input_options = CategorizedInputOptions(
        structure=structure,
        task=task,
        xc=xc,
        **options
    )
    
    # Generate VASP input files
    vif = VaspInputFiles(
        input_options, 
        overridden_incar_settings=user_incar_settings or {}
    )
    
    return VaspInputSet(
        incar=vif.incar,
        kpoints=vif.kpoints,
        poscar=vif.poscar,
        potcar=vif.potcar,
        structure=vif._structure,
        initial_structure=vif.initial_structure
    )


def get_available_tasks() -> list:
    """
    Get list of available task types.
    
    Returns:
        List of task name strings
    """
    return Task.name_list()


def get_available_xc() -> list:
    """
    Get list of available exchange-correlation functionals.
    
    Returns:
        List of XC functional name strings
    """
    return Xc.name_list()
