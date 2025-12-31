# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
# This file is originally developed by Naoki Tsunoda.
"""VESTA file generation module.

This module provides classes and utilities for generating VESTA
(Visualization for Electronic and STructural Analysis) files from
pymatgen Structure objects. VESTA files (.vesta) are used to visualize
crystal structures, isosurfaces, and other crystallographic data.

Key Components:
    - VestaFile: Main class that aggregates all VESTA blocks
    - Block classes (Title, Cellp, Struc, etc.): Individual VESTA file sections
    - add_density: Utility function to add volumetric data to existing VESTA files
    - calc_isurfs: Utility function to calculate isosurface values
"""
from __future__ import annotations

import itertools
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
from pymatgen.core import DummySpecies
from pymatgen.core.structure import Structure

from vise.analyzer.vesta.element_colors import atom_color
from vise.util.logger import get_logger
from vise.util.unit_conversion import au_to_angstrom

logger = get_logger(__name__)

# Type aliases for improved readability
RGBColor = Tuple[int, int, int]
"""RGB color tuple with values 0-255."""

Vector3D = Iterable[float]
"""3D vector as an iterable of floats."""


def replace_dummy_to_xx(target_str: str) -> str:
    """Replace pymatgen dummy species notation with VESTA notation.

    Pymatgen uses 'X0+' for dummy species, while VESTA uses 'XX'.

    Args:
        target_str: String containing pymatgen dummy species notation.

    Returns:
        String with 'X0+' replaced by 'XX'.
    """
    return target_str.replace("X0+", "XX", -1)


def val_to_str_line(values: Union[list, tuple]) -> str:
    """Convert a list or tuple of numerical values to a space-separated string.

    Each value is formatted to 6 decimal places.

    Args:
        values: Sequence of numerical values.

    Returns:
        Space-separated string of formatted values.

    Example:
        >>> val_to_str_line([1.0, 2.0, 3.0])
        '1.000000 2.000000 3.000000'
    """
    format_str = "{:.6f}"
    return " ".join([format_str.format(v) for v in values])


class VestaFile:
    """Main class for generating VESTA visualization files.

    This class aggregates all VESTA file blocks (Title, Cellp, Struc, etc.)
    and provides methods to generate the complete file content.

    Attributes:
        blocks: List of VESTA block objects that compose the file.
    """

    def __init__(
        self,
        structure: Structure,
        title: Optional[str] = None,
        vectors: Optional[Dict[int, Vector3D]] = None,
        vector_colors: Optional[List[RGBColor]] = None,
        bond_radius: float = 0.12,
        boundary: Optional[Iterable[float]] = None,
        show_label: bool = True,
    ) -> None:
        """Initialize a VestaFile with the given structure and options.

        If a name is set to the Site in structure (via properties), they are
        used for atom labels in the visualization.

        Args:
            structure: pymatgen Structure object to visualize.
            title: Custom title for the VESTA file. Defaults to chemical formula.
            vectors: Dictionary mapping site indices to vector directions.
                Example: {1: [1.0, 0.0, 0.0], 2: [0.0, 1.0, 0.0]}
            vector_colors: List of RGB color tuples for each vector.
                Values should be in range 0-255.
            bond_radius: Radius for bond visualization in Angstrom.
            boundary: Boundary limits as (a_min, a_max, b_min, b_max, c_min, c_max).
                Center of origin is determined from the mean of boundary.
                Example: (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5) means (0,0,0) as center.
            show_label: Whether to show atom labels in the visualization.
        """
        self.blocks = [
            Title(structure, title),
            Cellp(structure),
            Struc(structure),
            Bound(boundary),
            SBond(structure, bond_radius=bond_radius),
            SiteT(structure, show_label=show_label),
            # DummyAtomt(structure),  # Currently unused
            Vect(vectors, size=0.2, colors=vector_colors),
            Style(bond_radius=bond_radius),
        ]

    def __repr__(self) -> str:
        """Return the complete VESTA file content as a string."""
        block_representations = [repr(block) for block in self.blocks]
        return "\n\n".join(block_representations)

    def write_file(self, filename: Union[Path, str]) -> None:
        """Write the VESTA file to disk.

        Args:
            filename: Output filename. The '.vesta' extension is added
                automatically if not present.
        """
        filename_str = str(filename)
        if ".vesta" not in filename_str:
            filename_str = f"{filename_str}.vesta"
        with open(filename_str, "w") as vesta_file:
            vesta_file.write(repr(self))
        logger.info(f"{filename_str} is created.")


class Title:
    """TITLE block in VESTA files showing the visualization title.

    Attributes:
        vesta_version: VESTA format version string.
        header: Block header identifier.
        title: The title string to display.
    """

    vesta_version: str = "#VESTA_FORMAT_VERSION 3.5.0"
    header: str = "TITLE"

    def __init__(self, structure: Structure, title: Optional[str] = None) -> None:
        """Initialize the Title block.

        Args:
            structure: pymatgen Structure object.
            title: Custom title. Defaults to the structure's chemical formula.
        """
        self.title = title or structure.formula

    def __repr__(self) -> str:
        """Return the formatted TITLE block."""
        return "\n".join([self.vesta_version, self.header, self.title])


class Cellp:
    """CELLP block in VESTA files describing lattice constants.

    This block contains the six lattice parameters: a, b, c (lengths in Angstrom)
    and alpha, beta, gamma (angles in degrees).

    Example output:
        CELLP
        3.992878   3.992878   3.992878  90.000000  90.000000  90.000000

    Attributes:
        header: Block header identifier.
        cell_param: Tuple of (a, b, c, alpha, beta, gamma).
    """

    header: str = "CELLP"

    def __init__(self, structure: Structure) -> None:
        """Initialize the Cellp block from a structure.

        Args:
            structure: pymatgen Structure object.
        """
        self.cell_param = structure.lattice.lengths + structure.lattice.angles

    def __repr__(self) -> str:
        """Return the formatted CELLP block."""
        return "\n".join([self.header, val_to_str_line(self.cell_param)])


class Struc:
    """STRUC block in VESTA files showing atom sites.

    Each atom entry contains: index, element symbol, label, occupation,
    and fractional coordinates. Site indices are 1-based.

    Example output:
        STRUC
        1       Ba        Ba1               1.0000         0.500  0.500  0.500
        0 0 0 0 0

    Attributes:
        header: Block header identifier.
        separator: Block terminator line.
        zero_coord: Zero coordinate padding line for each site.
        str_coords: Formatted string containing all site coordinates.
    """

    header: str = "STRUC"
    separator: str = " 0 0 0 0 0 "
    zero_coord: str = " 0.0 0.0 0.0 "

    def __init__(self, structure: Structure) -> None:
        """Initialize the Struc block from a structure.

        Args:
            structure: pymatgen Structure object.
        """
        coord_lines: List[str] = []
        for site_idx, site in enumerate(structure, start=1):
            # Use site name from properties if available, otherwise generate one
            site_name = site.properties.get("name", None) or f"{site.species_string}{site_idx}"
            site_info = f"{site_idx} {site.species_string} {site_name} {1.0} "
            frac_coords_str = val_to_str_line(site.frac_coords)
            coord_lines.append(site_info + frac_coords_str)
            coord_lines.append(self.zero_coord)

        # Replace pymatgen dummy species notation "X0+" with VESTA notation "XX"
        self.str_coords = replace_dummy_to_xx("\n".join(coord_lines))

    def __repr__(self) -> str:
        """Return the formatted STRUC block."""
        return "\n".join([self.header, self.str_coords, self.separator])


class ImportDensity:
    """IMPORT_DENSITY block in VESTA files for volumetric data.

    This block specifies a volumetric data file (e.g., CHGCAR, PARCHG)
    to be imported and visualized.

    Attributes:
        smooth_param: Smoothing parameter for density visualization.
        header: Block header with smoothing parameter.
        prefix: Density scale prefix.
        string: Formatted volumetric file reference.
    """

    smooth_param: float = 0.1
    header: str = f"IMPORT_DENSITY {smooth_param}"
    prefix: str = "+1.000000"

    def __init__(self, volumetric_filename: str) -> None:
        """Initialize the ImportDensity block.

        Args:
            volumetric_filename: Name of the volumetric data file
                (e.g., 'CHGCAR', 'PARCHG').
        """
        self.string = f"{self.prefix} {volumetric_filename}"

    def __repr__(self) -> str:
        """Return the formatted IMPORT_DENSITY block."""
        return "\n".join([self.header, self.string])


class Bound:
    """BOUND block controlling the visualization boundary range.

    Defines the fractional coordinate limits for structure display.

    Attributes:
        header: Block header identifier.
        separator: Block terminator line.
        boundary: Formatted boundary string or None if not set.
    """

    header: str = "BOUND"
    separator: str = " 0 0 0 0 0 "

    def __init__(
        self,
        boundary: Optional[Union[tuple, list]] = (0, 1, 0, 1, 0, 1),
    ) -> None:
        """Initialize the Bound block.

        Args:
            boundary: Boundary limits as (a_min, a_max, b_min, b_max, c_min, c_max).
                Must contain exactly 6 values.

        Raises:
            ValueError: If boundary does not contain exactly 6 values.
        """
        if boundary:
            if len(boundary) != 6:
                raise ValueError("length of boundary must be 6")
            self.boundary = val_to_str_line(boundary)
        else:
            self.boundary = None

    def __repr__(self) -> str:
        """Return the formatted BOUND block, or empty string if no boundary."""
        if self.boundary:
            return "\n".join([self.header, self.boundary, self.separator])
        return ""


class SBond:
    """SBOND block controlling bond visualization between elements.

    Defines which atom pairs should be connected with bonds and their
    visualization parameters.

    Attributes:
        header: Block header identifier.
        separator: Block terminator line.
        bonds: Formatted string containing all bond definitions.

    Notes:
        boundary_mode = 0: "Do not search the atoms beyond the boundary"
    """

    header: str = "SBOND"
    separator: str = " 0 0 0 0 "

    def __init__(
        self,
        structure: Structure,
        bond_factor: float = 1.2,
        boundary_mode: int = 0,
        bond_radius: Optional[float] = None,
    ) -> None:
        """Initialize the SBond block from a structure.

        Args:
            structure: pymatgen Structure object.
            bond_factor: Factor to multiply ionic radii sum for bond length cutoff.
            boundary_mode: VESTA boundary mode (0 = do not search beyond boundary).
            bond_radius: Optional fixed bond radius for visualization.
        """
        bond_lines: List[str] = []
        element_pairs = itertools.permutations(structure.types_of_species, 2)

        for bond_idx, (elem1, elem2) in enumerate(element_pairs, start=1):
            try:
                bond_length = float(
                    elem1.average_ionic_radius + elem2.average_ionic_radius
                ) * bond_factor
            except AttributeError:
                # Skip pairs where ionic radius is not available (e.g., DummySpecies)
                continue

            bond_line = (
                f"{bond_idx} {elem1.symbol} {elem2.symbol} 0.0 {bond_length:5.4}  "
                f"0  {boundary_mode}  1  0  1"
            )
            if bond_radius:
                bond_line += f" {bond_radius}   2.000 161  33 246"
            bond_lines.append(bond_line)

        self.bonds = "\n".join(bond_lines)

    def __repr__(self) -> str:
        """Return the formatted SBOND block."""
        return "\n".join([self.header, self.bonds, self.separator])


class SiteT:
    """SITET block controlling atom site visualization properties.

    Defines properties like color, radius, and label visibility for each
    atomic site in the structure.

    Attributes:
        header: Block header identifier.
        separator: Block terminator line.
        sites: List of formatted site property strings.

    Notes:
        boundary_mode = 0: "Do not search the atoms beyond the boundary"
    """

    header: str = "SITET"
    separator: str = " 0 0 0 0 "

    def __init__(self, structure: Structure, show_label: bool = True) -> None:
        """Initialize the SiteT block from a structure.

        Args:
            structure: pymatgen Structure object.
            show_label: Whether to display atom labels.
        """
        self.sites: List[str] = []

        for site_idx, site in enumerate(structure, start=1):
            site_name = site.properties.get("name", None) or f"{site.species_string}{site_idx}"

            # Use gray color for dummy species, element-specific color otherwise
            if isinstance(site.specie, DummySpecies):
                rgb_str = "30 30 30"
            else:
                rgb_str = atom_color(site.species_string)

            # Reduce radius for "center" marker sites
            radius = 0.3 if site_name == "center" else 0.5

            # Determine label visibility
            if show_label:
                label_flag = 0 if site_name == "center" else 1
            else:
                label_flag = 0

            self.sites.append(f" {site_idx}  {site_name}  {radius} {rgb_str} {rgb_str} 204 {label_flag}")

    def __repr__(self) -> str:
        """Return the formatted SITET block."""
        return "\n".join([self.header] + self.sites + [self.separator])


class Vect:
    """VECTR/VECTT blocks for vector visualization.

    Displays directional arrows on specified atomic sites, useful for
    showing forces, displacements, or magnetic moments.

    Vector type options:
        - type 0: No penetration, no add_atomic_radius
        - type 1: Penetration, no add_atomic_radius
        - type 2: No penetration, add_atomic_radius
        - type 3: Penetration, add_atomic_radius

    Attributes:
        header_1: VECTR block header identifier.
        separator: Block terminator line.
        header_2: VECTT block header identifier.
        type: Vector visualization type.
        vectors: Formatted VECTR block content or None.
        vector_types: Formatted VECTT block content.
    """

    header_1: str = "VECTR"
    separator: str = " 0 0 0 0 0 "
    header_2: str = "VECTT"
    type: int = 2

    def __init__(
        self,
        vectors: Optional[Dict[int, Vector3D]],
        size: float = 0.5,
        colors: Optional[List[RGBColor]] = None,
    ) -> None:
        """Initialize the Vect block with vector definitions.

        Args:
            vectors: Dictionary mapping site indices to vector components.
                Example: {1: [1.0, 0.0, 0.0], 2: [0.0, 1.0, 0.0]}
            size: Vector arrow size.
            colors: List of RGB color tuples (0-255) for each vector.
                (0, 0, 0) = black. Defaults to all black if not specified.
        """
        if vectors:
            vec_lines: List[str] = []
            type_lines: List[str] = []

            # Default to black color for all vectors
            if colors is None:
                colors = [(0, 0, 0)] * len(vectors)

            for vec_idx, (site_idx, vec) in enumerate(vectors.items(), start=1):
                # VECTR block: vector index, components, site index, options
                vec_lines.append(f"{vec_idx} {val_to_str_line(vec)}")
                vec_lines.append(f"{site_idx}  0 0 0 0")
                vec_lines.append(self.separator)

                # VECTT block: vector index, size, RGB color, type
                color = colors[vec_idx - 1]
                color_str = f"{color[0]} {color[1]} {color[2]}"
                type_lines.append(f"{vec_idx} {size} {color_str} {self.type}")

            vec_lines.append(self.separator)
            type_lines.append(self.separator)

            self.vectors = "\n".join(vec_lines)
            self.vector_types = "\n".join(type_lines)
        else:
            self.vectors = None

    def __repr__(self) -> str:
        """Return the formatted VECTR/VECTT blocks, or empty string if no vectors."""
        if self.vectors:
            return "\n".join([
                self.header_1,
                self.vectors,
                "",
                self.header_2,
                self.vector_types,
            ])
        return ""


class DummyAtomt:
    """ATOMT block for controlling dummy atom visualization style.

    This class is currently not used in VestaFile but is retained for
    potential future use.

    Attributes:
        header: Block header identifier.
        str: Formatted atom style string or empty.
    """

    header: str = "ATOMT"

    def __init__(self, structure: Structure) -> None:
        """Initialize the DummyAtomt block.

        Args:
            structure: pymatgen Structure object.
        """
        if "X" in structure.symbol_set:
            self.str = "  1         Xx  0.2000  76  76  76  76  76  76 504"
        else:
            self.str = ""

    def __repr__(self) -> str:
        """Return the formatted ATOMT block, or empty string if no dummy atoms."""
        if self.str:
            return "\n".join([self.header, self.str])
        return ""


class Style:
    """STYLE block controlling overall visualization style.

    This block contains various visualization settings like vector amplitude,
    section parameters, atom display mode, and bond properties.

    Attributes:
        header: Block header identifier.
        Various style-related class attributes for formatting.
    """

    header: str = "STYLE"
    amp_prefix: str = "VECTS "
    atoms_prefix: str = "ATOMS "
    bondp_prefix: str = "BONDP"
    sects_prefix: str = "SECTS"
    sectp_prefix: str = "SECTP"
    ucolp_prefix: str = "UCOLP"

    all_bold_cell: str = "0  2  1.000   0   0   0"
    isurf: str = "1   1    12.0991 255 255   0 127 255"
    sect_param: int = 64
    sects: str = f" {sect_param}  1"
    sectp: str = "  1 0 0  0.00000E+00  0.00000E+00  -1.00000E+01  1.00000E+01"

    def __init__(self, bond_radius: float, is_ionic: bool = True) -> None:
        """Initialize the Style block.

        Args:
            bond_radius: Bond visualization radius in Angstrom.
            is_ionic: If True, use ionic display mode. Defaults to True.
        """
        self.amplitude: float = 1.0
        self.atoms = "1  0  1" if is_ionic else "0  0  1"
        self.bond_radius = bond_radius

    def __repr__(self) -> str:
        """Return the formatted STYLE block."""
        vct = f"{self.amp_prefix} {self.amplitude}"
        sects = f"{self.sects_prefix} {self.sects}"
        sectp = f"{self.sectp_prefix} \n {self.sectp}"
        ucol = f"{self.ucolp_prefix} \n {self.all_bold_cell}"
        atoms = f"{self.atoms_prefix} {self.atoms}"
        bondp = f"{self.bondp_prefix} \n   1  16  {self.bond_radius}  1.000 127 127 127"
        return "\n".join([self.header, vct, sects, sectp, ucol, atoms, bondp])


def add_density(
    original_vesta_text: str,
    isurfs: List[float],
    volumetric_filename: str,
) -> str:
    """Add volumetric density data to an existing VESTA file.

    This function modifies an existing VESTA file text to include volumetric
    data (e.g., charge density, partial charge) and isosurface definitions.

    Note:
        When 'CHG' is included in the volumetric file name, VESTA
        automatically divides the quantity by volume in cubic atomic units.

    Args:
        original_vesta_text: The original VESTA file content as a string.
        isurfs: List of isosurface values to display.
        volumetric_filename: Name of the volumetric data file (e.g., 'CHGCAR').

    Returns:
        Modified VESTA file content with density visualization added.
    """
    lines = original_vesta_text.split("\n")

    # Find insertion point after TITLE block
    title_idx = lines.index("TITLE")
    import_density_line = repr(ImportDensity(volumetric_filename))
    lines.insert(title_idx + 3, import_density_line + "\n")

    # Add ISURF block for isosurface definitions
    lines.append("\nISURF")
    for isurf_value in isurfs:
        lines.append(f"  1   1  {isurf_value}  0  0  255  50  50")
    lines.append("  0   0   0   0")

    return "\n".join(lines)


def calc_isurfs(
    border_fractions: List[float],
    is_chg: bool,
    volume: float,
) -> List[float]:
    """Calculate ISURF values for use in VESTA files.

    This function converts border fraction values to the appropriate
    isosurface values used by VESTA. It handles the unit conversion
    required when dealing with charge density files.

    Note:
        This calculation is valid only for light-weighted volumetric data.

    Args:
        border_fractions: List of fraction values defining isosurface levels.
        is_chg: If True, apply charge density unit conversion (Angstrom to a.u.).
        volume: Volume of the unit cell in cubic Angstroms.

    Returns:
        List of calculated ISURF values rounded to 5 decimal places.
    """
    # Scale fractions by the number of levels (max value = len(border_fractions))
    isurfs = np.array(border_fractions) * len(border_fractions)

    if is_chg:
        # VESTA uses atomic units for length in charge density visualization
        volume_in_au = volume / (au_to_angstrom**3)
        isurfs /= volume_in_au

    return np.round(isurfs, 5).tolist()
