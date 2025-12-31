# -*- coding: utf-8 -*-
#  Copyright (c) 2023. Distributed under the terms of the MIT License.
"""Coordinate sorting utilities for vise.

This module provides functions for sorting 3D coordinates based on
angular relationships, useful for crystal structure analysis.
"""

import numpy as np
from numpy.linalg import det


def sort_coords(coords: np.ndarray) -> np.ndarray:
    """Sort 3D coordinates by angle relative to the first coordinate.

    Sorts coordinates based on the angle between the ray from the center
    to each coordinate and the ray from the center to the first coordinate.
    The sorting is performed in the plane defined by the first two
    non-parallel vectors.

    Args:
        coords: Array of 3D coordinates with shape (n, 3) where n >= 2.
                Format: np.array([[x1, y1, z1], [x2, y2, z2], ...])

    Returns:
        Array of sorted coordinates with the same shape as input.

    Raises:
        ValueError: If coordinates are not 3-dimensional.

    Examples:
        >>> coords = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]])
        >>> sorted_coords = sort_coords(coords)
    """
    if len(coords[0]) != 3:
        raise ValueError("Only valid for 3D vectors")

    # Calculate center point
    center = np.average(coords, axis=0)
    relative_coords = coords - center

    # Find normal vector to the plane
    cross_product = np.cross(relative_coords[0], relative_coords[1])

    # Handle parallel vectors by using the third coordinate
    if abs(np.linalg.norm(cross_product)) < 1e-8 and len(relative_coords) > 2:
        cross_product = np.cross(relative_coords[0], relative_coords[2])

    plane_normal = cross_product / np.linalg.norm(cross_product)

    # Reference vector (normalized first coordinate)
    ref_vector = relative_coords[0] / np.linalg.norm(relative_coords[0])

    def angle_from_reference(index: int) -> float:
        """Calculate angle between reference vector and vector at index.

        Args:
            index: Index of the coordinate to compare.

        Returns:
            Angle in radians between the reference vector and the vector
            from center to coords[index].
        """
        v = relative_coords[index] / np.linalg.norm(relative_coords[index])
        matrix = np.concatenate(([ref_vector], [v], [plane_normal]), axis=0)
        determinant = det(matrix)
        dot_product = np.clip(np.dot(ref_vector, v), -1.0, 1.0)
        return np.arctan2(dot_product, determinant)

    # Sort indices by angle
    indices = list(range(len(coords)))
    indices.sort(key=angle_from_reference)

    return coords[indices]
