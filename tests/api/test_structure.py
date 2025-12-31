# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Tests for vise.api.structure module."""

import pytest
from pymatgen.core import Structure, Lattice

from vise.api.structure import (
    get_symmetry_info,
    get_primitive,
    get_conventional,
    SymmetryInfo,
)


@pytest.fixture
def cubic_structure():
    """Create a simple cubic NaCl structure."""
    lattice = Lattice.cubic(5.64)
    species = ["Na", "Na", "Na", "Na", "Cl", "Cl", "Cl", "Cl"]
    coords = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.5, 0.5],
    ]
    return Structure(lattice, species, coords)


@pytest.fixture
def simple_cubic():
    """Create a simple cubic structure (Pm-3m)."""
    lattice = Lattice.cubic(4.0)
    return Structure(lattice, ["Si"], [[0, 0, 0]])


class TestGetSymmetryInfo:
    """Tests for get_symmetry_info function."""

    def test_returns_symmetry_info(self, cubic_structure):
        """Test that function returns SymmetryInfo dataclass."""
        result = get_symmetry_info(cubic_structure)
        assert isinstance(result, SymmetryInfo)

    def test_space_group_symbol(self, cubic_structure):
        """Test that space group symbol is returned."""
        result = get_symmetry_info(cubic_structure)
        assert isinstance(result.space_group_symbol, str)
        assert len(result.space_group_symbol) > 0

    def test_space_group_number(self, cubic_structure):
        """Test that space group number is in valid range."""
        result = get_symmetry_info(cubic_structure)
        assert 1 <= result.space_group_number <= 230

    def test_nacl_is_fcc(self, cubic_structure):
        """Test that NaCl is identified as Fm-3m."""
        result = get_symmetry_info(cubic_structure)
        assert result.space_group_number == 225
        assert "Fm-3m" in result.space_group_symbol

    def test_point_group(self, cubic_structure):
        """Test that point group is returned."""
        result = get_symmetry_info(cubic_structure)
        assert isinstance(result.point_group, str)

    def test_crystal_system(self, cubic_structure):
        """Test that crystal system is returned."""
        result = get_symmetry_info(cubic_structure)
        assert result.crystal_system in [
            "cubic", "tetragonal", "orthorhombic", "hexagonal",
            "trigonal", "monoclinic", "triclinic"
        ]

    def test_volume(self, cubic_structure):
        """Test that volume is correctly calculated."""
        result = get_symmetry_info(cubic_structure)
        expected_volume = 5.64 ** 3
        assert abs(result.volume - expected_volume) < 0.01

    def test_lattice_abc(self, cubic_structure):
        """Test that lattice parameters are returned."""
        result = get_symmetry_info(cubic_structure)
        assert len(result.lattice_abc) == 3
        for param in result.lattice_abc:
            assert param > 0

    def test_lattice_angles(self, cubic_structure):
        """Test that lattice angles are returned."""
        result = get_symmetry_info(cubic_structure)
        assert len(result.lattice_angles) == 3
        # Cubic structure should have 90 degree angles
        for angle in result.lattice_angles:
            assert abs(angle - 90.0) < 0.1

    def test_custom_symprec(self, cubic_structure):
        """Test with custom symprec parameter."""
        result = get_symmetry_info(cubic_structure, symprec=0.1)
        assert isinstance(result, SymmetryInfo)

    def test_custom_angle_tolerance(self, cubic_structure):
        """Test with custom angle_tolerance parameter."""
        result = get_symmetry_info(cubic_structure, angle_tolerance=10.0)
        assert isinstance(result, SymmetryInfo)

    def test_simple_cubic(self, simple_cubic):
        """Test simple cubic structure."""
        result = get_symmetry_info(simple_cubic)
        assert result.space_group_number == 221  # Pm-3m


class TestGetPrimitive:
    """Tests for get_primitive function."""

    def test_returns_structure(self, cubic_structure):
        """Test that function returns a Structure."""
        result = get_primitive(cubic_structure)
        assert isinstance(result, Structure)

    def test_primitive_has_fewer_or_equal_atoms(self, cubic_structure):
        """Test that primitive cell has fewer or equal atoms."""
        result = get_primitive(cubic_structure)
        assert len(result) <= len(cubic_structure)

    def test_nacl_primitive_has_2_atoms(self, cubic_structure):
        """Test that NaCl primitive cell has 2 atoms."""
        result = get_primitive(cubic_structure)
        assert len(result) == 2

    def test_custom_symprec(self, cubic_structure):
        """Test with custom symprec."""
        result = get_primitive(cubic_structure, symprec=0.1)
        assert isinstance(result, Structure)


class TestGetConventional:
    """Tests for get_conventional function."""

    def test_returns_structure(self, cubic_structure):
        """Test that function returns a Structure."""
        result = get_conventional(cubic_structure)
        assert isinstance(result, Structure)

    def test_conventional_has_more_or_equal_atoms(self, cubic_structure):
        """Test that conventional cell has more or equal atoms than primitive."""
        primitive = get_primitive(cubic_structure)
        conventional = get_conventional(cubic_structure)
        assert len(conventional) >= len(primitive)

    def test_custom_symprec(self, cubic_structure):
        """Test with custom symprec."""
        result = get_conventional(cubic_structure, symprec=0.1)
        assert isinstance(result, Structure)


class TestSymmetryInfoDataclass:
    """Tests for SymmetryInfo dataclass."""

    def test_all_fields_present(self):
        """Test that all expected fields are present."""
        info = SymmetryInfo(
            space_group_symbol="Fm-3m",
            space_group_number=225,
            point_group="m-3m",
            crystal_system="cubic",
            bravais_lattice="cF",
            centering="F",
            volume=179.4,
            lattice_abc=(5.64, 5.64, 5.64),
            lattice_angles=(90.0, 90.0, 90.0),
            is_primitive=False
        )
        
        assert info.space_group_symbol == "Fm-3m"
        assert info.space_group_number == 225
        assert info.point_group == "m-3m"
        assert info.crystal_system == "cubic"
        assert info.bravais_lattice == "cF"
        assert info.centering == "F"
        assert info.volume == 179.4
        assert info.lattice_abc == (5.64, 5.64, 5.64)
        assert info.lattice_angles == (90.0, 90.0, 90.0)
        assert info.is_primitive is False
