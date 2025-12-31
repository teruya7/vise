Tutorial -- Input Set
---------------------

This tutorial shows how to generate VASP input files.

.. contents:: Table of Contents
   :local:
   :depth: 2

Input Files for Cell Relaxation
===============================

Let's assume we have a POSCAR file and are about to prepare INCAR, POTCAR, and KPOINTS files.
The sub-command for this is ``vasp_set`` (alias: ``vs``).

Important options are ``--task`` (``-t``) and ``--xc`` (``-x``) 
that determine the target task and exchange-correlation functional.
The defaults are structure optimization and PBE functional.

Generate input files:

.. code-block:: bash

    vise vs

Override POTCAR and INCAR settings:

.. code-block:: bash

    vise vs --potcar Mg_pv -uis ALGO All

You can also control settings via ``vise.yaml`` file, but command line arguments take priority.
See :doc:`tutorial_vise_yaml` for details.

Tasks
=====

List of INCAR settings for each task:

.. csv-table:: Tasks
   :file: task.csv
   :header-rows: 1

**Notes:**

- \*\*: Not set, uses VASP default
- \*\*1: Band count = valence electrons / 2 + unoccupied bands per element
- \*\*2: ISMEAR=-4 or -5 for insulators with ≥4 k-points, else ISMEAR=0
- \*\*3: IBRION=8 and NPAR≥2 combination is prohibited
- \*\*4: K-point multiplication factor (e.g., DOS doubles k-points)
- \*\*5: Band calculations require primitive cells from seekpath
- \*\*6: See `kpoints_mode.py <https://github.com/kumagai-group/vise/blob/master/vise/input_set/kpoints_mode.py>`_

**Common settings:**

- SIGMA=0.1, NELM=100, LASPH=True
- For DOS/dielectric: EMIN, EMAX, NEDOS are set based on band edges
- Gamma-centered mesh is forced when ISMEAR=-5

XC Functionals
==============

INCAR settings for each XC functional:

.. csv-table:: XC functional
   :file: xc.csv
   :header-rows: 1

**Note:** LDAUU and LDAUL parameters are in `u_parameter_set.yaml <https://github.com/kumagai-group/vise/blob/master/vise/input_set/datasets/u_parameter_set.yaml>`_

POTCAR Files
============

Default POTCAR selection is in `potcar_set.yaml <https://github.com/kumagai-group/vise/blob/master/vise/input_set/datasets/potcar_set.yaml>`_.

KPOINTS Files
=============

K-point mesh is proportional to reciprocal lattice constants.

For a cubic lattice with a=10Å:

- Reciprocal length: 2π/10
- With density 2.5Å: k-points = ceil(2π/10 × 2.5) = ceil(π/2) = 2

For body-centered orthorhombic/tetragonal systems, k-point counts are equalized 
using the geometric mean of reciprocal lattice constants.

Command Options
===============

--prev_dir
~~~~~~~~~~

Parse previous VASP calculations:

.. code-block:: bash

    vise vs --prev_dir ../relax

Extracts structure, charge, band edges, and magnetization.
Use with ``--file_transfer`` to copy/move/link files.

--options
~~~~~~~~~

Set generator options:

**IncarSettingsGenerator options:**

- ``charge``: System charge (default: 0.0)
- ``band_gap``: Band gap value
- ``vbm_cbm``: Band edge positions
- ``exchange_ratio``: Hybrid functional mixing (default: 0.25)
- ``set_hubbard_u``: Enable DFT+U
- ``auto_npar_kpar``: Auto-set parallelization (default: True)
- ``cutoff_energy``: Custom ENCUT
- ``is_magnetization``: Spin-polarized calculation

**StructureKpointsGenerator options:**

- ``kpt_density``: K-point density in Å
- ``gamma_centered``: Force Gamma centering
- ``only_even_num_kpts``: Use even k-point counts
- ``num_kpt_factor``: K-point multiplication factor

Example:

.. code-block:: bash

    vise vs --options cutoff_energy 1000 only_even_num_kpts True

See Also
========

- :doc:`tutorial_properties` - Analyze calculation results
- :doc:`cli_reference` - Full CLI reference
- :doc:`workflows` - Complete workflow examples
