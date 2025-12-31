# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Tests for vise.api.materials_project module."""

import pytest
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path

from vise.api.materials_project import (
    search_materials,
    get_structure_by_id,
    get_material_info,
    get_most_stable_structure,
    MaterialEntry,
    MaterialInfo
)


class TestSearchMaterials:
    """Tests for search_materials function."""

    def test_search_returns_list(self):
        """Test that search returns list of MaterialEntry."""
        with patch('mp_api.client.MPRester') as MockMPRester:
            mock_entry = Mock()
            mock_entry.material_id = "mp-149"
            mock_entry.energy_above_hull = 0.0
            mock_entry.symmetry = Mock(symbol="Fd-3m")
            mock_entry.band_gap = 1.1
            
            mock_rester = MagicMock()
            mock_rester.summary.search.return_value = [mock_entry]
            MockMPRester.return_value.__enter__.return_value = mock_rester
            
            result = search_materials("Si")
            
            assert isinstance(result, list)
            assert len(result) == 1
            assert isinstance(result[0], MaterialEntry)

    def test_search_result_attributes(self):
        """Test MaterialEntry attributes from search."""
        with patch('mp_api.client.MPRester') as MockMPRester:
            mock_entry = Mock()
            mock_entry.material_id = "mp-149"
            mock_entry.energy_above_hull = 0.0
            mock_entry.symmetry = Mock(symbol="Fd-3m")
            mock_entry.band_gap = 1.1
            
            mock_rester = MagicMock()
            mock_rester.summary.search.return_value = [mock_entry]
            MockMPRester.return_value.__enter__.return_value = mock_rester
            
            result = search_materials("Si")
            
            assert result[0].material_id == "mp-149"
            assert result[0].energy_above_hull == 0.0
            assert result[0].space_group == "Fd-3m"
            assert result[0].band_gap == 1.1

    def test_search_sorted_by_energy(self):
        """Test that results are sorted by energy above hull."""
        with patch('mp_api.client.MPRester') as MockMPRester:
            entry1 = Mock()
            entry1.material_id = "mp-1"
            entry1.energy_above_hull = 0.1
            entry1.symmetry = Mock(symbol="P1")
            entry1.band_gap = 1.0
            
            entry2 = Mock()
            entry2.material_id = "mp-2"
            entry2.energy_above_hull = 0.0  # Most stable
            entry2.symmetry = Mock(symbol="P2")
            entry2.band_gap = 1.5
            
            mock_rester = MagicMock()
            mock_rester.summary.search.return_value = [entry1, entry2]
            MockMPRester.return_value.__enter__.return_value = mock_rester
            
            result = search_materials("test")
            
            # Should be sorted, mp-2 first
            assert result[0].material_id == "mp-2"
            assert result[1].material_id == "mp-1"

    def test_search_empty_result(self):
        """Test empty search result."""
        with patch('mp_api.client.MPRester') as MockMPRester:
            mock_rester = MagicMock()
            mock_rester.summary.search.return_value = []
            MockMPRester.return_value.__enter__.return_value = mock_rester
            
            result = search_materials("NonexistentFormula")
            
            assert result == []


class TestGetStructureById:
    """Tests for get_structure_by_id function."""

    def test_returns_structure(self):
        """Test that function returns Structure."""
        with patch('mp_api.client.MPRester') as MockMPRester:
            from pymatgen.core import Structure, Lattice
            
            mock_structure = Structure(
                Lattice.cubic(5.43),
                ["Si", "Si"],
                [[0, 0, 0], [0.25, 0.25, 0.25]]
            )
            
            mock_rester = Mock()
            mock_rester.get_structure_by_material_id.return_value = mock_structure
            MockMPRester.return_value = mock_rester
            
            result = get_structure_by_id("mp-149")
            
            assert isinstance(result, Structure)

    def test_save_poscar(self, tmp_path):
        """Test save_poscar option."""
        with patch('mp_api.client.MPRester') as MockMPRester:
            from pymatgen.core import Structure, Lattice
            
            mock_structure = Structure(
                Lattice.cubic(5.43),
                ["Si", "Si"],
                [[0, 0, 0], [0.25, 0.25, 0.25]]
            )
            
            mock_rester = Mock()
            mock_rester.get_structure_by_material_id.return_value = mock_structure
            MockMPRester.return_value = mock_rester
            
            poscar_path = tmp_path / "POSCAR"
            result = get_structure_by_id(
                "mp-149", 
                save_poscar=True, 
                poscar_filename=str(poscar_path)
            )
            
            assert poscar_path.exists()


class TestGetMaterialInfo:
    """Tests for get_material_info function."""

    def test_returns_material_info(self):
        """Test that function returns MaterialInfo."""
        with patch('mp_api.client.MPRester') as MockMPRester:
            from pymatgen.core import Structure, Lattice
            
            mock_structure = Structure(
                Lattice.cubic(5.43),
                ["Si", "Si"],
                [[0, 0, 0], [0.25, 0.25, 0.25]]
            )
            
            mock_data = Mock()
            mock_data.total_magnetization = 0.0
            mock_data.band_gap = 1.1
            
            mock_rester = Mock()
            mock_rester.get_structure_by_material_id.return_value = mock_structure
            mock_rester.summary = Mock()
            mock_rester.summary.search.return_value = [mock_data]
            MockMPRester.return_value = mock_rester
            
            result = get_material_info("mp-149")
            
            assert isinstance(result, MaterialInfo)
            assert result.material_id == "mp-149"
            assert result.band_gap == 1.1


class TestGetMostStableStructure:
    """Tests for get_most_stable_structure function."""

    def test_returns_most_stable_structure(self):
        """Test that function returns most stable structure."""
        with patch('vise.api.materials_project.search_materials') as MockSearch, \
             patch('vise.api.materials_project.get_structure_by_id') as MockGetStruct:
            
            from pymatgen.core import Structure, Lattice
            
            mock_entry = MaterialEntry(
                material_id="mp-149",
                energy_above_hull=0.0,
                space_group="Fd-3m",
                band_gap=1.1
            )
            MockSearch.return_value = [mock_entry]
            
            mock_structure = Structure(
                Lattice.cubic(5.43),
                ["Si", "Si"],
                [[0, 0, 0], [0.25, 0.25, 0.25]]
            )
            MockGetStruct.return_value = mock_structure
            
            result = get_most_stable_structure("Si")
            
            assert isinstance(result, Structure)
            MockGetStruct.assert_called_once_with(
                "mp-149", save_poscar=False, poscar_filename="POSCAR"
            )

    def test_no_entries_returns_none(self):
        """Test that empty search returns None."""
        with patch('vise.api.materials_project.search_materials') as MockSearch:
            MockSearch.return_value = []
            
            result = get_most_stable_structure("NonexistentFormula")
            
            assert result is None


class TestMaterialEntry:
    """Tests for MaterialEntry class."""

    def test_material_entry_attributes(self):
        """Test MaterialEntry has expected attributes."""
        entry = MaterialEntry(
            material_id="mp-149",
            energy_above_hull=0.0,
            space_group="Fd-3m",
            band_gap=1.1
        )
        
        assert entry.material_id == "mp-149"
        assert entry.energy_above_hull == 0.0
        assert entry.space_group == "Fd-3m"
        assert entry.band_gap == 1.1


class TestMaterialInfo:
    """Tests for MaterialInfo class."""

    def test_material_info_attributes(self):
        """Test MaterialInfo has expected attributes."""
        from pymatgen.core import Structure, Lattice
        
        mock_structure = Structure(
            Lattice.cubic(5.43),
            ["Si", "Si"],
            [[0, 0, 0], [0.25, 0.25, 0.25]]
        )
        
        info = MaterialInfo(
            structure=mock_structure,
            material_id="mp-149",
            total_magnetization=0.0,
            band_gap=1.1
        )
        
        assert info.material_id == "mp-149"
        assert info.total_magnetization == 0.0
        assert info.band_gap == 1.1
        assert isinstance(info.structure, Structure)
