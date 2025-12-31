# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""Tests for vise.api.vasp_inputs module."""

import pytest
from pathlib import Path
from pymatgen.core import Structure, Lattice

from vise.api.vasp_inputs import (
    create_vasp_set,
    VaspInputSet,
    get_available_tasks,
    get_available_xc,
)
from vise.input_set.task import Task
from vise.input_set.xc import Xc


@pytest.fixture
def simple_structure():
    """Create a simple cubic Si structure."""
    lattice = Lattice.cubic(5.43)
    return Structure(lattice, ["Si", "Si"], [[0, 0, 0], [0.25, 0.25, 0.25]])


@pytest.fixture
def mgo_structure():
    """Create a simple MgO structure."""
    lattice = Lattice.cubic(4.21)
    return Structure(
        lattice, 
        ["Mg", "O"], 
        [[0, 0, 0], [0.5, 0.5, 0.5]]
    )





def potcar_available():
    """Check if POTCAR files are available."""
    try:
        from pymatgen.io.vasp import Potcar
        # Try to create a simple potcar
        Potcar(["Si"])
        return True
    except Exception:
        return False


class TestCreateVaspSet:
    """Tests for create_vasp_set function.
    
    Note: These tests require POTCAR files to be configured.
    """

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_returns_vasp_input_set(self, simple_structure):
        """Test that function returns VaspInputSet."""
        result = create_vasp_set(simple_structure)
        assert isinstance(result, VaspInputSet)

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_default_task_is_structure_opt(self, simple_structure):
        """Test that default task is structure_opt."""
        result = create_vasp_set(simple_structure)
        # IBRION should be set for structure optimization
        assert "IBRION" in result.incar

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_default_xc_is_pbe(self, simple_structure):
        """Test that default XC is PBE."""
        result = create_vasp_set(simple_structure)
        # PBE doesn't explicitly set GGA in INCAR (it's default)
        assert result.incar.get("GGA", "PE") in ["PE", None]

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_task_as_string(self, simple_structure):
        """Test task parameter as string."""
        result = create_vasp_set(simple_structure, task="structure_opt")
        assert isinstance(result, VaspInputSet)

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_task_as_enum(self, simple_structure):
        """Test task parameter as Task enum."""
        result = create_vasp_set(simple_structure, task=Task.structure_opt)
        assert isinstance(result, VaspInputSet)

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_xc_as_string(self, simple_structure):
        """Test xc parameter as string."""
        result = create_vasp_set(simple_structure, xc="pbe")
        assert isinstance(result, VaspInputSet)

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_xc_as_enum(self, simple_structure):
        """Test xc parameter as Xc enum."""
        result = create_vasp_set(simple_structure, xc=Xc.pbe)
        assert isinstance(result, VaspInputSet)

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_kpt_density_parameter(self, simple_structure):
        """Test kpt_density parameter."""
        result = create_vasp_set(simple_structure, kpt_density=2.0)
        assert isinstance(result, VaspInputSet)
        assert result.kpoints is not None

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_user_incar_settings(self, simple_structure):
        """Test user_incar_settings parameter."""
        result = create_vasp_set(
            simple_structure, 
            user_incar_settings={"EDIFF": 1e-7, "NELM": 200}
        )
        assert result.incar["EDIFF"] == 1e-7
        assert result.incar["NELM"] == 200

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_has_incar(self, simple_structure):
        """Test that result has INCAR."""
        result = create_vasp_set(simple_structure)
        assert result.incar is not None

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_has_kpoints(self, simple_structure):
        """Test that result has KPOINTS."""
        result = create_vasp_set(simple_structure)
        assert result.kpoints is not None

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_has_poscar(self, simple_structure):
        """Test that result has POSCAR."""
        result = create_vasp_set(simple_structure)
        assert result.poscar is not None

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_has_structure(self, simple_structure):
        """Test that result has structure."""
        result = create_vasp_set(simple_structure)
        assert result.structure is not None
        assert isinstance(result.structure, Structure)

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_has_initial_structure(self, simple_structure):
        """Test that result has initial_structure."""
        result = create_vasp_set(simple_structure)
        assert result.initial_structure is not None


class TestVaspInputSet:
    """Tests for VaspInputSet class."""

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_as_dict(self, simple_structure):
        """Test as_dict method."""
        inputs = create_vasp_set(simple_structure)
        d = inputs.as_dict()
        
        assert "incar" in d
        assert "kpoints" in d
        assert "structure" in d
        assert "num_atoms" in d

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_write_creates_directory(self, simple_structure, tmp_path):
        """Test that write creates directory if needed."""
        inputs = create_vasp_set(simple_structure)
        output_dir = tmp_path / "new_dir"
        
        inputs.write(output_dir)
        
        assert output_dir.exists()

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_write_creates_incar(self, simple_structure, tmp_path):
        """Test that write creates INCAR file."""
        inputs = create_vasp_set(simple_structure)
        inputs.write(tmp_path)
        
        assert (tmp_path / "INCAR").exists()

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_write_creates_kpoints(self, simple_structure, tmp_path):
        """Test that write creates KPOINTS file."""
        inputs = create_vasp_set(simple_structure)
        inputs.write(tmp_path)
        
        assert (tmp_path / "KPOINTS").exists()

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_write_creates_poscar(self, simple_structure, tmp_path):
        """Test that write creates POSCAR file."""
        inputs = create_vasp_set(simple_structure)
        inputs.write(tmp_path)
        
        assert (tmp_path / "POSCAR").exists()

    @pytest.mark.skipif(not potcar_available(), reason="POTCAR not available")
    def test_write_accepts_string_path(self, simple_structure, tmp_path):
        """Test that write accepts string path."""
        inputs = create_vasp_set(simple_structure)
        inputs.write(str(tmp_path))
        
        assert (tmp_path / "INCAR").exists()


class TestGetAvailableTasks:
    """Tests for get_available_tasks function."""

    def test_returns_list(self):
        """Test that function returns a list."""
        result = get_available_tasks()
        assert isinstance(result, list)

    def test_contains_structure_opt(self):
        """Test that list contains structure_opt."""
        result = get_available_tasks()
        # name_list() returns enum objects, check by converting
        result_names = [str(t) for t in result]
        assert "structure_opt" in result_names

    def test_contains_band(self):
        """Test that list contains band."""
        result = get_available_tasks()
        result_names = [str(t) for t in result]
        assert "band" in result_names

    def test_contains_dos(self):
        """Test that list contains dos."""
        result = get_available_tasks()
        result_names = [str(t) for t in result]
        assert "dos" in result_names


class TestGetAvailableXc:
    """Tests for get_available_xc function."""

    def test_returns_list(self):
        """Test that function returns a list."""
        result = get_available_xc()
        assert isinstance(result, list)

    def test_contains_pbe(self):
        """Test that list contains pbe."""
        result = get_available_xc()
        result_names = [str(x) for x in result]
        assert "pbe" in result_names

    def test_contains_hse(self):
        """Test that list contains hse."""
        result = get_available_xc()
        result_names = [str(x) for x in result]
        assert "hse" in result_names
