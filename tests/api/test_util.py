# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Tests for vise.api.util module - Unit tests for dataclasses and basic logic."""

import pytest
from unittest.mock import Mock
from pathlib import Path

from pymatgen.core import Structure, Lattice


class TestPhonopySetup:
    """Tests for PhonopySetup dataclass."""

    def test_phonopy_setup_attributes(self):
        """Test PhonopySetup has expected attributes."""
        from vise.api.util import PhonopySetup
        
        mock_supercell = Mock(spec=Structure)
        
        setup = PhonopySetup(
            supercell=mock_supercell,
            phonopy_input_file=Path("phonopy_input.json")
        )
        
        assert setup.supercell == mock_supercell
        assert setup.phonopy_input_file == Path("phonopy_input.json")

    def test_phonopy_setup_with_path_object(self):
        """Test PhonopySetup with Path object."""
        from vise.api.util import PhonopySetup
        
        lattice = Lattice.cubic(5.0)
        supercell = Structure(lattice, ["Si"] * 8, [
            [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
            [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], 
            [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]
        ])
        
        setup = PhonopySetup(
            supercell=supercell,
            phonopy_input_file=Path("/tmp/phonopy_input.json")
        )
        
        assert len(setup.supercell) == 8
        assert setup.phonopy_input_file.name == "phonopy_input.json"


class TestUtilFunctionImports:
    """Test that util functions are importable."""

    def test_import_make_atom_poscars(self):
        """Test make_atom_poscars is importable."""
        from vise.api.util import make_atom_poscars
        assert callable(make_atom_poscars)

    def test_import_make_spin_decomposed(self):
        """Test make_spin_decomposed_volumetric_files is importable."""
        from vise.api.util import make_spin_decomposed_volumetric_files
        assert callable(make_spin_decomposed_volumetric_files)

    def test_import_make_light_weight(self):
        """Test make_light_weight_volumetric_data is importable."""
        from vise.api.util import make_light_weight_volumetric_data
        assert callable(make_light_weight_volumetric_data)

    def test_import_create_vesta_file(self):
        """Test create_vesta_file is importable."""
        from vise.api.util import create_vesta_file
        assert callable(create_vesta_file)

    def test_import_create_phonon_setup(self):
        """Test create_phonon_setup is importable."""
        from vise.api.util import create_phonon_setup
        assert callable(create_phonon_setup)

    def test_import_analyze_phonon(self):
        """Test analyze_phonon is importable."""
        from vise.api.util import analyze_phonon
        assert callable(analyze_phonon)


class TestUtilFromApi:
    """Test that util functions are available from vise.api."""

    def test_import_from_api(self):
        """Test importing util functions from vise.api."""
        from vise.api import (
            make_atom_poscars,
            make_spin_decomposed_volumetric_files,
            make_light_weight_volumetric_data,
            create_vesta_file,
            create_phonon_setup,
            analyze_phonon,
            PhonopySetup,
        )
        
        assert callable(make_atom_poscars)
        assert callable(create_phonon_setup)
