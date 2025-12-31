# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

import pytest

from vise.util.unit_conversion import au_to_angstrom


class TestAuToAngstrom:
    """Tests for the au_to_angstrom constant."""

    def test_au_to_angstrom_value(self):
        """Test that au_to_angstrom has the expected value."""
        # The Bohr radius is approximately 0.529177 Angstroms
        assert 0.5291 < au_to_angstrom < 0.5292

    def test_au_to_angstrom_precise_value(self):
        """Test au_to_angstrom against known physical constant."""
        # More precise check - Bohr radius ≈ 0.52917721067 Å
        assert abs(au_to_angstrom - 0.52917721067) < 1e-8

    def test_au_to_angstrom_is_positive(self):
        """Test that the conversion factor is positive."""
        assert au_to_angstrom > 0

    def test_conversion_applied(self):
        """Test applying the conversion factor."""
        # 1 atomic unit = 0.529... Angstrom
        one_au_in_angstrom = 1.0 * au_to_angstrom
        assert abs(one_au_in_angstrom - 0.5291772) < 1e-6

        # 10 atomic units
        ten_au_in_angstrom = 10.0 * au_to_angstrom
        assert abs(ten_au_in_angstrom - 5.291772) < 1e-5
