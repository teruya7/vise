# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

import pytest

from vise.input_set.kpoints import ViseKpoints


class TestViseKpoints:
    """Tests for the ViseKpoints class."""

    def test_str_single_kpoint_automatic(self):
        """Test __str__ for automatic single k-point mesh."""
        # This test would require a real structure
        # ViseKpoints inherits from Kpoints, so we test via from_dict instead
        pass

    def test_from_dict(self):
        """Test creating ViseKpoints from a dictionary."""
        d = {
            "comment": "Test KPOINTS",
            "nkpoints": 0,
            "generation_style": "Gamma",
            "kpoints": [[4, 4, 4]],
            "usershift": [0.0, 0.0, 0.0],
            "kpts_weights": None,
            "coord_type": None,
            "labels": None,
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        assert kpoints.comment == "Test KPOINTS"
        assert kpoints.num_kpts == 0
        assert kpoints.style.name == "Gamma"

    def test_str_gamma_centered(self):
        """Test __str__ for Gamma-centered automatic mesh."""
        d = {
            "comment": "Gamma-centered mesh",
            "nkpoints": 0,
            "generation_style": "Gamma",
            "kpoints": [[4, 4, 4]],
            "usershift": [0.0, 0.0, 0.0],
            "kpts_weights": None,
            "coord_type": None,
            "labels": None,
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        result = str(kpoints)
        
        assert "Gamma-centered mesh" in result
        assert "Gamma" in result
        assert "4 4 4" in result

    def test_str_monkhorst_pack(self):
        """Test __str__ for Monkhorst-Pack mesh."""
        d = {
            "comment": "Monkhorst-Pack mesh",
            "nkpoints": 0,
            "generation_style": "Monkhorst",
            "kpoints": [[6, 6, 6]],
            "usershift": [0.0, 0.0, 0.0],
            "kpts_weights": None,
            "coord_type": None,
            "labels": None,
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        result = str(kpoints)
        
        assert "Monkhorst-Pack mesh" in result
        assert "Monkhorst" in result
        assert "6 6 6" in result

    def test_str_with_shift(self):
        """Test __str__ for mesh with non-zero shift."""
        d = {
            "comment": "Mesh with shift",
            "nkpoints": 0,
            "generation_style": "Gamma",
            "kpoints": [[4, 4, 4]],
            "usershift": [0.5, 0.5, 0.5],
            "kpts_weights": None,
            "coord_type": None,
            "labels": None,
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        result = str(kpoints)
        
        # Non-zero shift should be printed
        assert "0.5 0.5 0.5" in result

    def test_str_explicit_kpoints(self):
        """Test __str__ for explicit k-points list."""
        d = {
            "comment": "Explicit k-points",
            "nkpoints": 3,
            "generation_style": "Reciprocal",
            "kpoints": [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0]],
            "usershift": [0.0, 0.0, 0.0],
            "kpts_weights": [1, 2, 4],
            "coord_type": None,
            "labels": None,
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        result = str(kpoints)
        
        assert "Explicit k-points" in result
        assert "3" in result  # num_kpts

    def test_str_line_mode(self):
        """Test __str__ for line mode (band structure)."""
        d = {
            "comment": "Line mode for band structure",
            "nkpoints": 20,
            "generation_style": "Line_mode",
            "kpoints": [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], 
                        [0.5, 0.0, 0.0], [0.5, 0.5, 0.0]],
            "usershift": [0.0, 0.0, 0.0],
            "kpts_weights": None,
            "coord_type": "Reciprocal",
            "labels": ["\\Gamma", "X", "X", "M"],
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        result = str(kpoints)
        
        assert "Line mode" in result
        assert "Line_mode" in result
        # Labels should be included
        assert "\\Gamma" in result or "Gamma" in result

    def test_str_explicit_with_labels(self):
        """Test __str__ for explicit k-points with labels."""
        d = {
            "comment": "Explicit with labels",
            "nkpoints": 2,
            "generation_style": "Reciprocal",
            "kpoints": [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            "usershift": [0.0, 0.0, 0.0],
            "kpts_weights": [1, 8],
            "coord_type": None,
            "labels": ["Gamma", "R"],
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        result = str(kpoints)
        
        # Labels should appear in output
        assert "Gamma" in result
        assert "R" in result

    def test_inheritance_from_pymatgen_kpoints(self):
        """Test that ViseKpoints inherits from pymatgen Kpoints."""
        from pymatgen.io.vasp.sets import Kpoints
        
        d = {
            "comment": "Test",
            "nkpoints": 0,
            "generation_style": "Gamma",
            "kpoints": [[4, 4, 4]],
            "usershift": [0.0, 0.0, 0.0],
            "kpts_weights": None,
            "coord_type": None,
            "labels": None,
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        
        assert isinstance(kpoints, Kpoints)

    def test_as_dict(self):
        """Test that as_dict returns valid dictionary."""
        d = {
            "comment": "Test KPOINTS",
            "nkpoints": 0,
            "generation_style": "Gamma",
            "kpoints": [[4, 4, 4]],
            "usershift": [0.0, 0.0, 0.0],
            "kpts_weights": None,
            "coord_type": None,
            "labels": None,
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints"
        }
        kpoints = ViseKpoints.from_dict(d)
        result = kpoints.as_dict()
        
        assert "comment" in result
        assert result["comment"] == "Test KPOINTS"
