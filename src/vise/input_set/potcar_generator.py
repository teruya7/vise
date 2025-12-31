# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""POTCAR file generation for VASP.

This module provides utilities for generating POTCAR files with
appropriate pseudopotential selections for each element.
"""

from typing import Dict, List, Optional

from pymatgen.io.vasp.sets import Potcar

from vise.defaults import defaults
from vise.error import ViseError
from vise.input_set.datasets.potcar_set import PotcarSet
from vise.input_set.xc import Xc
from vise.util.logger import get_logger

logger = get_logger(__name__)


def generate_potcar(
    symbol_list: List[str],
    xc: Xc,
    potcar_set: PotcarSet = defaults.potcar_set,
    overridden_potcar: Optional[Dict[str, str]] = None,
) -> Potcar:
    """Generate a POTCAR file for the given elements.

    Args:
        symbol_list: List of element symbols (e.g., ["Si", "O"]).
        xc: Exchange-correlation functional (determines GGA vs LDA).
        potcar_set: POTCAR selection preset (default from vise defaults).
        overridden_potcar: Override specific elements' POTCAR choices.

    Returns:
        pymatgen Potcar object ready to write.

    Raises:
        ViseNoPotcarError: If POTCAR not found for an element.

    Examples:
        >>> potcar = generate_potcar(["Si", "O"], Xc.pbe)
        >>> potcar.write_file("POTCAR")
    """
    potcar_dict = potcar_set.overridden_potcar_dict(overridden_potcar)

    try:
        potcar_symbols = [potcar_dict[el] for el in symbol_list]
    except KeyError as e:
        raise ViseNoPotcarError(f"No POTCAR found for element: {e}")

    # Select appropriate functional library
    functional = defaults.lda_potcar if xc == Xc.lda else defaults.gga_potcar

    return Potcar(potcar_symbols, functional=functional)


class ViseNoPotcarError(ViseError):
    """Exception raised when POTCAR is not found for an element."""

    pass
