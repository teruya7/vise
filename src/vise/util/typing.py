# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
"""Type aliases for vise.

This module defines common type aliases used throughout the vise package
for better type hinting and code readability.
"""

from typing import List, Tuple, Union

# 3D coordinates as a tuple of floats (x, y, z)
Coords = Tuple[float, float, float]

# Generalized 3D coordinates that can be either a tuple or a list
GenCoords = Union[Tuple[float, float, float], List[float]]
