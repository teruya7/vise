====================
Preparing Structures
====================

In this tutorial, we show how to prepare POSCAR file via the Materials Project (MP) database.

.. contents:: Table of Contents
   :local:
   :depth: 2

Prerequisites
=============

Before using Materials Project integration, configure your API key.
See :doc:`installation` for details.

Quick verification:

.. code-block:: python

    from pymatgen.io.vasp.sets import VaspInputSet
    from pymatgen.core import Structure, Lattice
    from pymatgen.io.vasp.sets import MPRelaxSet
    s = Structure(Lattice.cubic(1), ["H", "He"], [[0.0]*3, [0.5]*3])
    vasp_set = MPRelaxSet(s)
    vasp_set.write_input(".")

If VASP files are not created, check your pymatgen configuration.

Getting POSCAR from Materials Project
=====================================

Command Line
~~~~~~~~~~~~

Use the :code:`get_poscar` (= :code:`gp`) sub-command:

.. code-block:: bash

    # By Materials Project ID
    vise gp -m mp-2857

    # By formula (gets most stable structure)
    vise gp -f SrTiO3

Python API
~~~~~~~~~~

.. code-block:: python

    from vise.api import materials_project

    # Search by formula
    entries = materials_project.search_materials("TiO2")
    for entry in entries[:5]:
        print(f"{entry.material_id}: {entry.space_group}")

    # Get structure by ID
    struct = materials_project.get_structure_by_id("mp-2857", save_poscar=True)

    # Get most stable structure
    struct = materials_project.get_most_stable_structure("SrTiO3")

Checking Structure Symmetry
===========================

After obtaining POSCAR, verify the structure:

.. code-block:: bash

    vise si

    # Output example:
    # Space group: Fm-3m (225)
    # Crystal system: cubic
    # Point group: m-3m

Get primitive or conventional cells:

.. code-block:: bash

    # Primitive cell
    vise si --primitive

    # Conventional cell
    vise si -c

Next Steps
==========

- :doc:`generating_inputs` - Generate VASP input files
- :doc:`cli_commands` - Full CLI command reference
