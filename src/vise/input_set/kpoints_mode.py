# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""K-point generation modes for VASP calculations.

This module defines the KpointsMode enum for specifying how k-points
should be generated for different calculation types.
"""

from vise.util.enum import ExtendedEnum


class KpointsMode(ExtendedEnum):
    """K-point mesh generation modes.

    Determines how k-points are generated for VASP calculations.

    Attributes:
        band: Generate k-points along high-symmetry band paths using
              seekpath. The primitive cell is used and returned with
              the band structure calculation setup.
        primitive: Uniform k-point mesh based on the standardized
                   primitive cell. Structure is converted to primitive
                   if necessary. Centering follows crystallographic conventions.
        uniform: Uniform k-point mesh on the given lattice. The mesh
                 is shifted along perpendicular directions for 90-degree
                 angles. Useful for supercell calculations.
    """

    band = "band"
    primitive = "primitive"
    uniform = "uniform"

    @property
    def band_or_primitive(self) -> bool:
        """Check if mode uses primitive cell symmetry analysis."""
        return self in (self.band, self.primitive)
