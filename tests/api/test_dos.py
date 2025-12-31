# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Tests for vise.api.dos module - Unit tests for dataclasses."""

import pytest
from unittest.mock import Mock

from vise.api.dos import DosAnalysisResult


class TestDosAnalysisResult:
    """Tests for DosAnalysisResult class."""

    def test_result_attributes(self):
        """Test result has expected attributes."""
        mock_dos_data = Mock()
        result = DosAnalysisResult(
            dos_data=mock_dos_data,
            band_gap=2.0,
            vbm=-0.5,
            cbm=1.5,
            is_metal=False,
            efermi=5.0
        )
        
        assert result.band_gap == 2.0
        assert result.vbm == -0.5
        assert result.cbm == 1.5
        assert result.is_metal is False
        assert result.efermi == 5.0
        assert result.dos_data == mock_dos_data

    def test_metallic_result(self):
        """Test metallic result attributes."""
        mock_dos_data = Mock()
        result = DosAnalysisResult(
            dos_data=mock_dos_data,
            band_gap=None,
            vbm=None,
            cbm=None,
            is_metal=True,
            efermi=5.0
        )
        
        assert result.is_metal is True
        assert result.band_gap is None
        assert result.efermi == 5.0

    def test_save_json(self, tmp_path):
        """Test save_json method."""
        mock_dos_data = Mock()
        mock_dos_data.to_json_file = Mock()
        
        result = DosAnalysisResult(
            dos_data=mock_dos_data,
            band_gap=2.0,
            vbm=-0.5,
            cbm=1.5,
            is_metal=False,
            efermi=5.0
        )
        
        output_file = tmp_path / "dos.json"
        result.save_json(str(output_file))
        
        mock_dos_data.to_json_file.assert_called_once_with(str(output_file))

    def test_semiconductor_cbm_vbm_relation(self):
        """Test CBM > VBM for semiconductor."""
        mock_dos_data = Mock()
        result = DosAnalysisResult(
            dos_data=mock_dos_data,
            band_gap=2.0,
            vbm=-0.5,
            cbm=1.5,
            is_metal=False,
            efermi=0.5
        )
        
        assert result.cbm > result.vbm
        assert abs(result.cbm - result.vbm - result.band_gap) < 0.001

    def test_default_json_filename(self):
        """Test default JSON filename."""
        mock_dos_data = Mock()
        mock_dos_data.to_json_file = Mock()
        
        result = DosAnalysisResult(
            dos_data=mock_dos_data,
            band_gap=2.0,
            vbm=-0.5,
            cbm=1.5,
            is_metal=False,
            efermi=5.0
        )
        
        result.save_json()
        mock_dos_data.to_json_file.assert_called_once()
