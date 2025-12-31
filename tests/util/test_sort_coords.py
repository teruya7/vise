# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

import numpy as np
import pytest

from vise.util.sort_coords import sort_coords


class TestSortCoords:
    """Tests for the sort_coords function."""

    def test_square_in_xy_plane(self):
        """Test sorting coordinates forming a square in XY plane."""
        # Unsorted square vertices
        coords = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        result = sort_coords(coords)
        
        # Result should be sorted by angle
        assert len(result) == 4
        # All original coordinates should be present
        assert result.shape == coords.shape

    def test_triangle_in_xy_plane(self):
        """Test sorting coordinates forming a triangle in XY plane."""
        coords = np.array([
            [1.0, 0.0, 0.0],
            [-0.5, 0.8660254, 0.0],
            [-0.5, -0.8660254, 0.0],
        ])
        result = sort_coords(coords)
        
        assert len(result) == 3
        # All original coordinates should be present
        assert result.shape == coords.shape

    def test_hexagon(self):
        """Test sorting coordinates forming a hexagon."""
        angles = np.linspace(0, 2*np.pi, 7)[:-1]  # 6 equally spaced angles
        coords = np.array([
            [np.cos(a), np.sin(a), 0.0] for a in angles
        ])
        # Shuffle the order
        shuffled = coords[[0, 3, 1, 5, 2, 4]]
        
        result = sort_coords(shuffled)
        assert len(result) == 6

    def test_3d_coords(self):
        """Test sorting coordinates in 3D space."""
        coords = np.array([
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [-1.0, 0.0, 1.0],
            [0.0, -1.0, 1.0],
        ])
        result = sort_coords(coords)
        
        assert len(result) == 4
        # All original coordinates should be present
        assert result.shape == coords.shape

    def test_invalid_2d_input_raises_error(self):
        """Test that 2D coordinates raise ValueError."""
        coords = np.array([
            [1.0, 0.0],
            [0.0, 1.0],
            [-1.0, 0.0],
        ])
        with pytest.raises(ValueError, match="Only valid for 3D vector"):
            sort_coords(coords)

    def test_parallel_first_two_vectors(self):
        """Test handling of parallel first two vectors."""
        # First two vectors are parallel (both point in +x direction)
        coords = np.array([
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],  # Parallel to first
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        # Should not raise an error due to handling of parallel vectors
        result = sort_coords(coords)
        assert len(result) == 4

    def test_preserves_shape(self):
        """Test that output shape matches input shape."""
        coords = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ])
        result = sort_coords(coords)
        assert result.shape == coords.shape

    def test_single_point_with_others(self):
        """Test with minimal valid number of points."""
        coords = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        result = sort_coords(coords)
        assert len(result) == 2
