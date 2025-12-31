# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Tests for vise.api.band module - Unit tests for dataclasses and basic logic."""

import pytest
from unittest.mock import Mock

from vise.api.band import BandAnalysisResult


class TestBandAnalysisResult:
    """Tests for BandAnalysisResult class."""

    def test_result_attributes(self):
        """Test result has expected attributes."""
        mock_plot_info = Mock()
        result = BandAnalysisResult(
            plot_info=mock_plot_info,
            band_gap=1.5,
            vbm=-0.5,
            cbm=1.0,
            is_metal=False
        )
        
        assert result.band_gap == 1.5
        assert result.vbm == -0.5
        assert result.cbm == 1.0
        assert result.is_metal is False
        assert result.plot_info == mock_plot_info

    def test_metallic_result(self):
        """Test metallic result attributes."""
        mock_plot_info = Mock()
        result = BandAnalysisResult(
            plot_info=mock_plot_info,
            band_gap=None,
            vbm=None,
            cbm=None,
            is_metal=True
        )
        
        assert result.is_metal is True
        assert result.band_gap is None
        assert result.vbm is None
        assert result.cbm is None

    def test_save_json(self, tmp_path):
        """Test save_json method."""
        mock_plot_info = Mock()
        mock_plot_info.to_json_file = Mock()
        
        result = BandAnalysisResult(
            plot_info=mock_plot_info,
            band_gap=1.5,
            vbm=-0.5,
            cbm=1.0,
            is_metal=False
        )
        
        output_file = tmp_path / "band.json"
        result.save_json(str(output_file))
        
        mock_plot_info.to_json_file.assert_called_once_with(str(output_file))

    def test_semiconductor_dataclass(self):
        """Test semiconductor with typical values."""
        mock_plot_info = Mock()
        result = BandAnalysisResult(
            plot_info=mock_plot_info,
            band_gap=2.5,
            vbm=4.2,
            cbm=6.7,
            is_metal=False
        )
        
        assert result.band_gap == 2.5
        assert result.cbm - result.vbm == 2.5

    def test_default_json_filename(self, tmp_path):
        """Test default JSON filename."""
        mock_plot_info = Mock()
        mock_plot_info.to_json_file = Mock()
        
        result = BandAnalysisResult(
            plot_info=mock_plot_info,
            band_gap=1.5,
            vbm=-0.5,
            cbm=1.0,
            is_metal=False
        )
        
        result.save_json()
        mock_plot_info.to_json_file.assert_called_once()
