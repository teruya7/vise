# -*- coding: utf-8 -*-
#  Copyright (c) 2023. Distributed under the terms of the MIT License.
"""Brillouin zone plotting data structures.

This module provides data classes for storing Brillouin zone visualization
information, including zone faces, high-symmetry points, and band paths.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from monty.json import MSONable

# Type aliases for clarity
CartesianCoord = List[float]  # [x, y, z] in Cartesian coordinates
FractionalCoord = List[float]  # [a, b, c] in fractional coordinates
Face = List[CartesianCoord]  # List of vertices forming a face
Path = List[CartesianCoord]  # List of points forming a path


@dataclass
class BZPlotInfo(MSONable):
    """Brillouin zone visualization data.

    Contains all information needed to render a 3D Brillouin zone plot
    with high-symmetry points and band paths.

    Attributes:
        faces: List of BZ faces, each as a list of vertex coordinates
               in Cartesian space.
        labels: High-symmetry point labels with coordinates.
               Format: {"Γ": {"cart": [0,0,0], "frac": [0,0,0]}, ...}
        band_paths: Optional band structure paths as connected points.
        rec_lat_vec: Reciprocal lattice vectors (3x3 matrix).

    Examples:
        >>> info = BZPlotInfo(faces=[...], labels={"X": {...}})
        >>> # Use with matplotlib or other 3D plotting libraries
    """

    faces: List[Face]
    labels: Dict[str, Dict[str, CartesianCoord]]
    band_paths: Optional[List[Path]] = None
    rec_lat_vec: Optional[List[CartesianCoord]] = None
