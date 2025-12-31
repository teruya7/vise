# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""VASP-specific band edge property extraction.

This module provides utilities for extracting band edge properties
from VASP calculation outputs.
"""

from collections import defaultdict
from typing import Dict, Tuple

import numpy as np
from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Outcar, Procar, Vasprun

from vise.analyzer.band_edge_properties import BandEdge, BandEdgeProperties


class VaspBandEdgeProperties(BandEdgeProperties):
    """Band edge analyzer for VASP calculations.

    Extracts VBM and CBM information from VASP vasprun.xml and OUTCAR.

    Attributes:
        structure: Final structure from VASP calculation.
    """

    def __init__(
        self,
        vasprun: Vasprun,
        outcar: Outcar,
        integer_criterion: float = 0.1,
    ) -> None:
        """Initialize from VASP output files.

        Args:
            vasprun: Parsed vasprun.xml.
            outcar: Parsed OUTCAR.
            integer_criterion: Threshold for integer occupation check.
        """
        total_mag = outcar.total_mag if outcar.total_mag else 0.0
        self.structure = vasprun.final_structure

        super().__init__(
            eigenvalues=eigenvalues_from_vasprun(vasprun),
            nelect=outcar.nelect,
            magnetization=total_mag,
            kpoint_coords=vasprun.actual_kpoints,
            integer_criterion=integer_criterion,
            is_non_collinear=vasprun.parameters.get("LNONCOLLINEAR", False),
        )


def edge_orbital_contributions(
    procar: Procar,
    structure: Structure,
    vbm_info: BandEdge,
    cbm_info: BandEdge,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]]]:
    """Calculate orbital contributions at band edges.

    Args:
        procar: Parsed PROCAR with orbital projections.
        structure: Crystal structure.
        vbm_info: VBM band edge information.
        cbm_info: CBM band edge information.

    Returns:
        Tuple of (vbm_composition, cbm_composition) dictionaries.
        Each maps element symbol to orbital contributions.

    Raises:
        ValueError: If orbital decomposition format is unrecognized.
    """
    elements = structure.composition.elements

    # Group atom indices by element
    element_indices = []
    for elem in elements:
        element_indices.append(
            [i for i, site in enumerate(structure) if site.specie == elem]
        )

    # Extract projections for VBM and CBM
    vbm_proj = procar.data[vbm_info.spin][vbm_info.kpoint_index][vbm_info.band_index]
    cbm_proj = procar.data[cbm_info.spin][cbm_info.kpoint_index][cbm_info.band_index]

    # Determine orbital groupings based on projection array size
    num_orbitals = len(vbm_proj[0])
    if num_orbitals == 3:
        orbital_indices = [[0], [1], [2]]
    elif num_orbitals == 4:
        orbital_indices = [[0], [1], [2], [3]]
    elif num_orbitals == 9:
        # s, p (3), d (5)
        orbital_indices = [[0], [1, 2, 3], [4, 5, 6, 7, 8]]
    elif num_orbitals == 16:
        # s, p (3), d (5), f (7)
        orbital_indices = [[0], [1, 2, 3], [4, 5, 6, 7, 8],
                           [9, 10, 11, 12, 13, 14, 15]]
    else:
        raise ValueError(f"Unrecognized orbital count: {num_orbitals}")

    vbm_comp: Dict[str, Dict[str, float]] = defaultdict(dict)
    cbm_comp: Dict[str, Dict[str, float]] = defaultdict(dict)

    orbital_names = ["s", "p", "d", "f"]
    for elem, elem_idxs in zip(elements, element_indices):
        for orb_name, orb_idxs in zip(orbital_names, orbital_indices):
            vbm_comp[str(elem)][orb_name] = vbm_proj[np.ix_(elem_idxs, orb_idxs)].sum()
            cbm_comp[str(elem)][orb_name] = cbm_proj[np.ix_(elem_idxs, orb_idxs)].sum()

    return dict(vbm_comp), dict(cbm_comp)


def eigenvalues_from_vasprun(vasprun: Vasprun) -> Dict[Spin, np.ndarray]:
    """Extract eigenvalues from vasprun, dropping occupations.

    Args:
        vasprun: Parsed vasprun.xml.

    Returns:
        Dictionary mapping Spin to eigenvalue array [kpoint][band].
    """
    return {spin: e[:, :, 0] for spin, e in vasprun.eigenvalues.items()}
