# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest

from vise.input_set.fft_grids import vasp_grid, RYTOEV, AUTOA, PI


class TestVaspGrid:
    """Tests for the vasp_grid function."""

    def test_vasp_grids_normal_precision(self):
        """Test FFT grid calculation with Normal precision."""
        actual = vasp_grid(encut=250, lattice_length=10.0, symprec="Normal")
        assert actual == 39

    def test_vasp_grids_high_precision(self):
        """Test FFT grid calculation with High precision (factor=4)."""
        actual = vasp_grid(encut=250, lattice_length=10.0, symprec="High")
        # High precision uses factor=4 instead of 3
        assert actual > 39  # Should be larger than normal

    def test_vasp_grids_accurate_precision(self):
        """Test FFT grid calculation with Accurate precision (factor=4)."""
        actual = vasp_grid(encut=250, lattice_length=10.0, symprec="Accurate")
        # Accurate precision uses factor=4
        assert actual > 39

    def test_vasp_grids_high_encut(self):
        """Test FFT grid calculation with high ENCUT."""
        actual = vasp_grid(encut=600, lattice_length=10.0, symprec="Normal")
        # Higher ENCUT should result in larger grid
        assert actual > 39

    def test_vasp_grids_small_lattice(self):
        """Test FFT grid calculation with small lattice constant."""
        actual = vasp_grid(encut=250, lattice_length=5.0, symprec="Normal")
        # Smaller lattice should result in smaller grid
        assert actual < 39

    def test_vasp_grids_large_lattice(self):
        """Test FFT grid calculation with large lattice constant."""
        actual = vasp_grid(encut=250, lattice_length=20.0, symprec="Normal")
        # Larger lattice should result in larger grid
        assert actual > 39

    def test_vasp_grids_case_insensitive(self):
        """Test that symprec is case insensitive for first character."""
        normal_lower = vasp_grid(encut=250, lattice_length=10.0, symprec="normal")
        normal_upper = vasp_grid(encut=250, lattice_length=10.0, symprec="Normal")
        assert normal_lower == normal_upper

    def test_vasp_grids_high_case_variants(self):
        """Test High precision with different case variants."""
        high1 = vasp_grid(encut=250, lattice_length=10.0, symprec="high")
        high2 = vasp_grid(encut=250, lattice_length=10.0, symprec="High")
        high3 = vasp_grid(encut=250, lattice_length=10.0, symprec="HIGH")
        assert high1 == high2 == high3

    def test_returns_integer(self):
        """Test that the function returns an integer."""
        result = vasp_grid(encut=300, lattice_length=8.0, symprec="Normal")
        assert isinstance(result, int)

    def test_positive_result(self):
        """Test that the result is always positive."""
        result = vasp_grid(encut=100, lattice_length=3.0, symprec="Normal")
        assert result > 0


class TestConstants:
    """Tests for the physical constants defined in fft_grids.py."""

    def test_rytoev_value(self):
        """Test RYTOEV constant (Rydberg to eV)."""
        # 1 Rydberg ≈ 13.6 eV
        assert 13.5 < RYTOEV < 13.7

    def test_autoa_value(self):
        """Test AUTOA constant (atomic unit to Angstrom)."""
        # 1 Bohr ≈ 0.529 Angstrom
        assert 0.52 < AUTOA < 0.54

    def test_pi_value(self):
        """Test PI constant."""
        import math
        assert abs(PI - math.pi) < 1e-10


