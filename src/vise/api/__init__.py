# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""
vise.api - Clean Python API for VASP input generation and analysis

This module provides a user-friendly Python API for vise functionality,
separate from the CLI layer. It allows programmatic access to:

- Structure symmetry analysis
- VASP input file generation
- Band structure analysis
- DOS analysis
- Dielectric function analysis
- Band edge properties and effective mass
- Materials Project integration

Usage:
    >>> from vise.api import structure, vasp_inputs, band, dos
    >>> from pymatgen.core import Structure
    >>>
    >>> # Structure analysis
    >>> struct = Structure.from_file("POSCAR")
    >>> info = structure.get_symmetry_info(struct)
    >>> print(f"Space group: {info.space_group_symbol}")
    >>>
    >>> # VASP input generation
    >>> inputs = vasp_inputs.create_vasp_set(struct, task="structure_opt")
    >>> inputs.write("./calculation")
    >>>
    >>> # Band structure analysis (requires VASP output)
    >>> result = band.analyze_band("vasprun.xml")
    >>> result.save_plot("band.pdf")
    >>>
    >>> # DOS analysis
    >>> dos_result = dos.analyze_dos("vasprun.xml", "OUTCAR")
    >>> dos_result.save_plot("dos.pdf")
"""

# Structure API
from vise.api.structure import (
    get_symmetry_info,
    get_primitive,
    get_conventional,
    SymmetryInfo,
)

# VASP inputs API
from vise.api.vasp_inputs import (
    create_vasp_set,
    VaspInputSet,
    get_available_tasks,
    get_available_xc,
)

# Band structure API
from vise.api.band import (
    analyze_band,
    BandAnalysisResult,
)

# DOS API
from vise.api.dos import (
    analyze_dos,
    DosAnalysisResult,
)

# Dielectric function API
from vise.api.dielectric import (
    analyze_dielectric,
    load_dielectric_from_csv,
    DielectricAnalysisResult,
)

# Band edge and effective mass API
from vise.api.band_edge import (
    get_band_edge_properties,
    calculate_effective_mass,
    BandEdgeResult,
    EffectiveMassResult,
    KpointInfo,
)

# Materials Project API
from vise.api.materials_project import (
    search_materials,
    get_structure_by_id,
    get_material_info,
    get_most_stable_structure,
    MaterialEntry,
    MaterialInfo,
)

# Utility API
from vise.api.util import (
    make_atom_poscars,
    make_spin_decomposed_volumetric_files,
    make_light_weight_volumetric_data,
    create_vesta_file,
    create_phonon_setup,
    analyze_phonon,
    PhonopySetup,
)

__all__ = [
    # Structure API
    "get_symmetry_info",
    "get_primitive",
    "get_conventional",
    "SymmetryInfo",
    # VASP inputs API
    "create_vasp_set",
    "VaspInputSet",
    "get_available_tasks",
    "get_available_xc",
    # Band structure API
    "analyze_band",
    "BandAnalysisResult",
    # DOS API
    "analyze_dos",
    "DosAnalysisResult",
    # Dielectric function API
    "analyze_dielectric",
    "load_dielectric_from_csv",
    "DielectricAnalysisResult",
    # Band edge and effective mass API
    "get_band_edge_properties",
    "calculate_effective_mass",
    "BandEdgeResult",
    "EffectiveMassResult",
    "KpointInfo",
    # Materials Project API
    "search_materials",
    "get_structure_by_id",
    "get_material_info",
    "get_most_stable_structure",
    "MaterialEntry",
    "MaterialInfo",
    # Utility API
    "make_atom_poscars",
    "make_spin_decomposed_volumetric_files",
    "make_light_weight_volumetric_data",
    "create_vesta_file",
    "create_phonon_setup",
    "analyze_phonon",
    "PhonopySetup",
]
