# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Tests for vise.api.dielectric module - Unit tests for dataclasses."""

import pytest
from unittest.mock import Mock

from vise.api.dielectric import DielectricAnalysisResult


class TestDielectricAnalysisResult:
    """Tests for DielectricAnalysisResult class."""

    def test_result_attributes(self):
        """Test result has expected attributes."""
        mock_data = Mock()
        result = DielectricAnalysisResult(diele_func_data=mock_data)
        
        assert result.diele_func_data == mock_data

    def test_save_json(self, tmp_path):
        """Test save_json method."""
        mock_data = Mock()
        mock_data.to_json_file = Mock()
        
        result = DielectricAnalysisResult(diele_func_data=mock_data)
        
        output_file = tmp_path / "diele.json"
        result.save_json(str(output_file))
        
        mock_data.to_json_file.assert_called_once_with(str(output_file))

    def test_save_csv(self, tmp_path):
        """Test save_csv method."""
        mock_data = Mock()
        mock_data.to_csv_file = Mock()
        
        result = DielectricAnalysisResult(diele_func_data=mock_data)
        
        output_file = tmp_path / "diele.csv"
        result.save_csv(str(output_file))
        
        mock_data.to_csv_file.assert_called_once_with(str(output_file))

    def test_default_json_filename(self):
        """Test default JSON filename."""
        mock_data = Mock()
        mock_data.to_json_file = Mock()
        
        result = DielectricAnalysisResult(diele_func_data=mock_data)
        result.save_json()
        
        mock_data.to_json_file.assert_called_once()

    def test_default_csv_filename(self):
        """Test default CSV filename."""
        mock_data = Mock()
        mock_data.to_csv_file = Mock()
        
        result = DielectricAnalysisResult(diele_func_data=mock_data)
        result.save_csv()
        
        mock_data.to_csv_file.assert_called_once()
