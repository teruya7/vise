# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Lattice centering types and transformation matrices.

This module defines the Centering enum representing crystallographic
lattice centering types (P, A, C, R, I, F) and provides transformation
matrices between conventional and primitive cells.
"""

import numpy as np

from vise.util.enum import ExtendedEnum


class Centering(ExtendedEnum):
    """Crystallographic lattice centering types.

    Each centering type provides methods to convert between conventional
    and primitive cell representations.

    Attributes:
        P: Primitive (no centering)
        A: A-face centered
        C: C-face centered
        R: Rhombohedral
        I: Body-centered (Innenzentriert)
        F: Face-centered (all faces)
    """

    P = "P"
    A = "A"
    C = "C"
    R = "R"
    I = "I"
    F = "F"

    @property
    def conv_to_primitive(self) -> np.ndarray:
        """Transformation matrix from conventional to primitive cell.

        Note:
            These matrices are transposed from those in spglib documentation
            due to different matrix conventions.
            See: https://spglib.github.io/spglib/definition.html

        Returns:
            3x3 numpy array representing the transformation matrix.

        Raises:
            NotImplementedError: If matrix is not defined for this centering.
        """
        if self is Centering.P:
            return np.eye(3)
        elif self is Centering.A:
            return np.array([
                [1,    0,    0],
                [0,  1/2,  1/2],
                [0, -1/2,  1/2]
            ])
        elif self is Centering.C:
            return np.array([
                [1/2, -1/2, 0],
                [1/2,  1/2, 0],
                [  0,    0, 1]
            ])
        elif self is Centering.R:
            return np.array([
                [ 2/3,  1/3,  1/3],
                [-1/3,  1/3,  1/3],
                [-1/3, -2/3,  1/3]
            ])
        elif self is Centering.I:
            return np.array([
                [-1/2,  1/2,  1/2],
                [ 1/2, -1/2,  1/2],
                [ 1/2,  1/2, -1/2]
            ])
        elif self is Centering.F:
            return np.array([
                [  0,  1/2, 1/2],
                [1/2,    0, 1/2],
                [1/2,  1/2,   0]
            ])
        else:
            raise NotImplementedError(f"Matrix not defined for {self}")

    @property
    def primitive_to_conv(self) -> np.ndarray:
        """Transformation matrix from primitive to conventional cell.

        Returns:
            3x3 numpy array (inverse of conv_to_primitive).
        """
        return np.linalg.inv(self.conv_to_primitive)

    @property
    def conv_multiplicity(self) -> int:
        """Number of primitive cells per conventional cell.

        Returns:
            Integer multiplicity (1 for P, 2 for A/C/I, 3 for R, 4 for F).
        """
        return int(round(1 / np.linalg.det(self.conv_to_primitive)))
