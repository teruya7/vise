# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Tests for vise.api.band_edge module - Unit tests for dataclasses."""

import pytest
from unittest.mock import Mock

from vise.api.band_edge import (
    BandEdgeResult,
    EffectiveMassResult,
    KpointInfo
)


class TestBandEdgeResult:
    """Tests for BandEdgeResult class."""

    def test_result_attributes(self):
        """Test BandEdgeResult has expected attributes."""
        vbm_info = KpointInfo(
            kpoint_index=0,
            kpoint_coords=(0.0, 0.0, 0.0),
            energy=-0.5
        )
        cbm_info = KpointInfo(
            kpoint_index=0,
            kpoint_coords=(0.0, 0.0, 0.0),
            energy=1.5
        )
        
        result = BandEdgeResult(
            band_gap=2.0,
            vbm_info=vbm_info,
            cbm_info=cbm_info,
            is_direct=True,
            is_metal=False,
            efermi=0.5
        )
        
        assert result.band_gap == 2.0
        assert result.is_direct is True
        assert result.is_metal is False
        assert result.efermi == 0.5
        assert result.vbm_info == vbm_info
        assert result.cbm_info == cbm_info

    def test_str_semiconductor_direct(self):
        """Test string representation for direct gap semiconductor."""
        vbm_info = KpointInfo(
            kpoint_index=0,
            kpoint_coords=(0.0, 0.0, 0.0),
            energy=-0.5
        )
        cbm_info = KpointInfo(
            kpoint_index=0,
            kpoint_coords=(0.0, 0.0, 0.0),
            energy=1.5
        )
        
        result = BandEdgeResult(
            band_gap=2.0,
            vbm_info=vbm_info,
            cbm_info=cbm_info,
            is_direct=True,
            is_metal=False,
            efermi=0.5
        )
        
        s = str(result)
        assert "2.0" in s
        assert "direct" in s

    def test_str_semiconductor_indirect(self):
        """Test string representation for indirect gap semiconductor."""
        vbm_info = KpointInfo(
            kpoint_index=0,
            kpoint_coords=(0.0, 0.0, 0.0),
            energy=-0.5
        )
        cbm_info = KpointInfo(
            kpoint_index=5,
            kpoint_coords=(0.5, 0.5, 0.5),
            energy=1.2
        )
        
        result = BandEdgeResult(
            band_gap=1.7,
            vbm_info=vbm_info,
            cbm_info=cbm_info,
            is_direct=False,
            is_metal=False,
            efermi=0.35
        )
        
        s = str(result)
        assert "1.7" in s
        assert "indirect" in s

    def test_str_metallic(self):
        """Test string representation for metal."""
        result = BandEdgeResult(
            band_gap=None,
            vbm_info=None,
            cbm_info=None,
            is_direct=False,
            is_metal=True,
            efermi=5.0
        )
        
        s = str(result)
        assert "Metallic" in s
        assert "5.0" in s

    def test_metallic_result_attributes(self):
        """Test metallic result attributes."""
        result = BandEdgeResult(
            band_gap=None,
            vbm_info=None,
            cbm_info=None,
            is_direct=False,
            is_metal=True,
            efermi=5.0
        )
        
        assert result.is_metal is True
        assert result.band_gap is None
        assert result.vbm_info is None


class TestKpointInfo:
    """Tests for KpointInfo class."""

    def test_kpoint_info_attributes(self):
        """Test KpointInfo has expected attributes."""
        info = KpointInfo(
            kpoint_index=5,
            kpoint_coords=(0.5, 0.5, 0.0),
            energy=-0.5
        )
        
        assert info.kpoint_index == 5
        assert info.kpoint_coords == (0.5, 0.5, 0.0)
        assert info.energy == -0.5

    def test_kpoint_info_gamma(self):
        """Test Gamma point info."""
        info = KpointInfo(
            kpoint_index=0,
            kpoint_coords=(0.0, 0.0, 0.0),
            energy=0.0
        )
        
        assert info.kpoint_index == 0
        assert info.kpoint_coords == (0.0, 0.0, 0.0)

    def test_kpoint_info_x_point(self):
        """Test X point info."""
        info = KpointInfo(
            kpoint_index=10,
            kpoint_coords=(0.5, 0.0, 0.5),
            energy=-1.2
        )
        
        assert info.kpoint_coords == (0.5, 0.0, 0.5)


class TestEffectiveMassResult:
    """Tests for EffectiveMassResult class."""

    def test_result_attributes(self):
        """Test EffectiveMassResult has expected attributes."""
        mock_eff_mass = Mock()
        result = EffectiveMassResult(
            effective_mass=mock_eff_mass,
            temperature=300.0,
            concentrations=[1e18, 1e19, 1e20]
        )
        
        assert result.temperature == 300.0
        assert result.concentrations == [1e18, 1e19, 1e20]
        assert result.effective_mass == mock_eff_mass

    def test_save_json(self, tmp_path):
        """Test save_json method."""
        mock_eff_mass = Mock()
        mock_eff_mass.to_json_file = Mock()
        
        result = EffectiveMassResult(
            effective_mass=mock_eff_mass,
            temperature=300.0,
            concentrations=[1e18]
        )
        
        output_file = tmp_path / "eff_mass.json"
        result.save_json(str(output_file))
        
        mock_eff_mass.to_json_file.assert_called_once_with(str(output_file))

    def test_str(self):
        """Test string representation."""
        mock_eff_mass = Mock()
        mock_eff_mass.__str__ = Mock(return_value="Effective mass: 0.5 m_e")
        
        result = EffectiveMassResult(
            effective_mass=mock_eff_mass,
            temperature=300.0,
            concentrations=[1e18]
        )
        
        s = str(result)
        assert "Effective mass" in s

    def test_default_json_filename(self):
        """Test default JSON filename."""
        mock_eff_mass = Mock()
        mock_eff_mass.to_json_file = Mock()
        
        result = EffectiveMassResult(
            effective_mass=mock_eff_mass,
            temperature=300.0,
            concentrations=[1e18]
        )
        
        result.save_json()
        mock_eff_mass.to_json_file.assert_called_once()
