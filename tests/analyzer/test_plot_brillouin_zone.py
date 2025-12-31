# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

import pytest

from vise.analyzer.plot_brillouin_zone import BZPlotInfo


class TestBZPlotInfo:
    """Tests for the BZPlotInfo dataclass."""

    def test_basic_initialization(self):
        """Test basic initialization with required fields."""
        faces = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]]
        labels = {"Gamma": {"cart": [0.0, 0.0, 0.0], "frac": [0.0, 0.0, 0.0]}}
        
        bz_info = BZPlotInfo(faces=faces, labels=labels)
        
        assert bz_info.faces == faces
        assert bz_info.labels == labels
        assert bz_info.band_paths is None
        assert bz_info.rec_lat_vec is None

    def test_full_initialization(self):
        """Test initialization with all fields."""
        faces = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]]
        labels = {
            "Gamma": {"cart": [0.0, 0.0, 0.0], "frac": [0.0, 0.0, 0.0]},
            "X": {"cart": [0.5, 0.0, 0.0], "frac": [0.7514, 0.0, 0.0]}
        }
        band_paths = [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]]
        rec_lat_vec = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        
        bz_info = BZPlotInfo(
            faces=faces, 
            labels=labels,
            band_paths=band_paths,
            rec_lat_vec=rec_lat_vec
        )
        
        assert bz_info.faces == faces
        assert bz_info.labels == labels
        assert bz_info.band_paths == band_paths
        assert bz_info.rec_lat_vec == rec_lat_vec

    def test_msonable_as_dict(self):
        """Test that BZPlotInfo can be serialized to dict via MSONable."""
        faces = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]]
        labels = {"Gamma": {"cart": [0.0, 0.0, 0.0], "frac": [0.0, 0.0, 0.0]}}
        
        bz_info = BZPlotInfo(faces=faces, labels=labels)
        d = bz_info.as_dict()
        
        assert "@module" in d
        assert "@class" in d
        assert d["@class"] == "BZPlotInfo"
        assert d["faces"] == faces
        assert d["labels"] == labels

    def test_msonable_from_dict(self):
        """Test that BZPlotInfo can be deserialized from dict via MSONable."""
        faces = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]]
        labels = {"X": {"cart": [0.5, 0.0, 0.0], "frac": [0.7514, 0.0, 0.0]}}
        
        original = BZPlotInfo(faces=faces, labels=labels)
        d = original.as_dict()
        restored = BZPlotInfo.from_dict(d)
        
        assert restored.faces == original.faces
        assert restored.labels == original.labels

    def test_msonable_roundtrip(self):
        """Test MSONable serialization round-trip."""
        faces = [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0]]
        ]
        labels = {
            "Gamma": {"cart": [0.0, 0.0, 0.0], "frac": [0.0, 0.0, 0.0]},
            "X": {"cart": [0.5, 0.0, 0.0], "frac": [0.7514, 0.0, 0.0]},
            "M": {"cart": [0.5, 0.5, 0.0], "frac": [0.5, 0.5, 0.0]}
        }
        band_paths = [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0]]]
        rec_lat_vec = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        
        original = BZPlotInfo(
            faces=faces, 
            labels=labels,
            band_paths=band_paths,
            rec_lat_vec=rec_lat_vec
        )
        d = original.as_dict()
        restored = BZPlotInfo.from_dict(d)
        
        assert restored.faces == original.faces
        assert restored.labels == original.labels
        assert restored.band_paths == original.band_paths
        assert restored.rec_lat_vec == original.rec_lat_vec

    def test_empty_faces(self):
        """Test with empty faces list."""
        bz_info = BZPlotInfo(faces=[], labels={})
        assert bz_info.faces == []
        assert bz_info.labels == {}

    def test_multiple_faces(self):
        """Test with multiple faces (typical BZ has many faces)."""
        # Simplified hexagonal BZ faces
        faces = [
            [[0.0, 0.0, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.5]],
            [[0.0, 0.0, -0.5], [0.5, 0.0, -0.5], [0.5, 0.5, -0.5]],
            [[0.5, 0.0, 0.5], [0.5, 0.0, -0.5], [0.5, 0.5, 0.0]],
        ]
        labels = {"K": {"cart": [0.333, 0.333, 0.0], "frac": [0.5, 0.5, 0.0]}}
        
        bz_info = BZPlotInfo(faces=faces, labels=labels)
        
        assert len(bz_info.faces) == 3

    def test_labels_structure(self):
        """Test that labels have correct structure."""
        labels = {
            "Gamma": {"cart": [0.0, 0.0, 0.0], "frac": [0.0, 0.0, 0.0]},
            "X": {"cart": [0.5, 0.0, 0.0], "frac": [0.7514, 0.0, 0.0]}
        }
        
        bz_info = BZPlotInfo(faces=[], labels=labels)
        
        assert "Gamma" in bz_info.labels
        assert "cart" in bz_info.labels["Gamma"]
        assert "frac" in bz_info.labels["Gamma"]
        assert len(bz_info.labels["Gamma"]["cart"]) == 3
