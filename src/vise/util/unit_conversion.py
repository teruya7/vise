# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
"""Unit conversion constants for vise.

This module provides physical constants used for unit conversions,
particularly between atomic units and Angstroms.
"""

from scipy.constants import physical_constants

# Conversion factor from atomic units (Bohr radius) to Angstroms
# 1 Bohr = 0.5291772109 Å
AU_TO_ANGSTROM: float = physical_constants["atomic unit of length"][0] * 1e10

# Backward compatibility alias
au_to_angstrom = AU_TO_ANGSTROM
