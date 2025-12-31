# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Experimental dielectric function data.

This module provides access to experimental dielectric function data
from literature for comparison with calculated values.
"""

from pathlib import Path
from typing import Dict

import pandas as pd

from vise.util.logger import get_logger

logger = get_logger(__name__)

# Directory containing experimental data CSV files
_DATA_DIR = Path(__file__).parent

# Experimental band gaps for reference materials (in eV)
BAND_GAPS: Dict[str, float] = {
    "GaAs": 1.41,
    "Si": 1.17,
}


class ExpDieleFunc:
    """Experimental dielectric function data reader.

    Loads experimental optical data from CSV files for comparison
    with calculated dielectric functions.

    Attributes:
        material: Material name (e.g., "Si", "GaAs").
        dataframe: Pandas DataFrame with optical data.

    Examples:
        >>> exp = ExpDieleFunc("Si")
        >>> exp.band_gap
        1.17
    """

    def __init__(self, material: str) -> None:
        """Initialize with material name.

        Args:
            material: Material identifier matching CSV filename.

        Raises:
            FileNotFoundError: If no CSV file exists for the material.
        """
        self.material = material
        try:
            self.dataframe = pd.read_csv(
                _DATA_DIR / f"{material}.csv", sep=r"\s+"
            )
        except FileNotFoundError:
            logger.warning(f"CSV file for {material} does not exist.")
            raise

    @property
    def energies(self) -> pd.Series:
        """Photon energies in eV."""
        return self.dataframe.energy

    @property
    def dielectric_real(self) -> pd.Series:
        """Real part of dielectric function."""
        return self.dataframe.real_part

    @property
    def dielectric_imag(self) -> pd.Series:
        """Imaginary part of dielectric function."""
        return self.dataframe.imaginary_part

    @property
    def absorption_coeff(self) -> pd.Series:
        """Absorption coefficient in cm⁻¹."""
        return self.dataframe.absorption

    @property
    def band_gap(self) -> float:
        """Experimental band gap in eV."""
        return BAND_GAPS[self.material]

    @property
    def reference(self) -> Dict[str, str]:
        """Literature reference for the data."""
        return {
            "title": "THE HANDBOOK ON OPTICAL CONSTANTS OF SEMICONDUCTORS",
            "author": "Sadao Adachi",
            "publisher": "World Scientific",
            "doi": "doi.org/10.1142/8480",
            "year": "2012",
        }


# Backward compatibility alias
band_gap = BAND_GAPS
