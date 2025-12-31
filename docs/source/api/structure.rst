==================
vise.api.structure
==================

.. module:: vise.api.structure
   :synopsis: Crystal structure symmetry analysis

This module provides functions for analyzing crystal structure symmetry
and transforming between primitive and conventional cells.

Quick Example
-------------

::

    from pymatgen.core import Structure
    from vise.api import structure

    # Load a structure
    struct = Structure.from_file("POSCAR")

    # Get symmetry information
    info = structure.get_symmetry_info(struct)
    print(f"Space group: {info.space_group_symbol} ({info.space_group_number})")
    print(f"Crystal system: {info.crystal_system}")
    print(f"Point group: {info.point_group}")
    print(f"Bravais lattice: {info.bravais_lattice}")
    print(f"Volume: {info.volume:.4f} Å³")

    # Get primitive cell
    primitive = structure.get_primitive(struct)
    primitive.to(filename="POSCAR_primitive")

    # Get conventional cell  
    conventional = structure.get_conventional(struct)
    conventional.to(filename="POSCAR_conventional")

Data Classes
------------

SymmetryInfo
~~~~~~~~~~~~

.. autoclass:: SymmetryInfo
   :members:
   :undoc-members:

Functions
---------

get_symmetry_info
~~~~~~~~~~~~~~~~~

.. autofunction:: get_symmetry_info

get_primitive
~~~~~~~~~~~~~

.. autofunction:: get_primitive

get_conventional
~~~~~~~~~~~~~~~~

.. autofunction:: get_conventional

See Also
--------

- :doc:`vasp_inputs` - Generate VASP input files from structures
- CLI command: ``vise vs`` for symmetry analysis
