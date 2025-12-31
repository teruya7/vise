# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""Exchange-correlation functional definitions.

This module defines the Xc enum representing different exchange-correlation
functionals supported in VASP calculations, from LDA through hybrid functionals.
"""

from monty.json import MSONable

from vise.util.enum import ExtendedEnum


class Xc(MSONable, ExtendedEnum):
    """Exchange-correlation functional types for VASP.

    Supports LDA, GGA, meta-GGA, and hybrid functionals commonly used
    in DFT calculations.

    Attributes:
        pbe: Perdew-Burke-Ernzerhof GGA functional
        pbesol: PBE revised for solids
        lda: Local Density Approximation
        scan: Strongly Constrained and Appropriately Normed meta-GGA
        r2scan: Regularized-restored SCAN
        pbe0: PBE-based hybrid with 25% exact exchange
        hse: Heyd-Scuseria-Ernzerhof screened hybrid
    """

    pbe = "pbe"
    pbesol = "pbesol"
    lda = "lda"
    scan = "scan"
    r2scan = "r2scan"
    pbe0 = "pbe0"
    hse = "hse"

    @classmethod
    def from_string(cls, name: str) -> "Xc":
        """Create Xc from string, handling alternative names.

        Args:
            name: Functional name (case-sensitive).

        Returns:
            Corresponding Xc enum member.

        Examples:
            >>> Xc.from_string("pbe")
            pbe
            >>> Xc.from_string("perdew-zunger81")
            lda
        """
        # Handle alternative LDA names
        if name == "perdew-zunger81":
            return cls.lda
        return super().from_string(name)

    @property
    def is_lda_or_gga(self) -> bool:
        """Check if functional is LDA or GGA type."""
        return self in (self.pbe, self.pbesol, self.lda)

    @property
    def is_metagga(self) -> bool:
        """Check if functional is a meta-GGA."""
        return self in (self.scan, self.r2scan)

    @property
    def is_hybrid_functional(self) -> bool:
        """Check if functional is a hybrid (includes exact exchange)."""
        return self in (self.pbe0, self.hse)

    @property
    def is_local_or_semilocal(self) -> bool:
        """Check if functional is local or semilocal (LDA, GGA, or meta-GGA)."""
        return self.is_lda_or_gga or self.is_metagga

    @property
    def is_nonlocal(self) -> bool:
        """Check if functional includes nonlocal exchange."""
        return self.is_hybrid_functional
