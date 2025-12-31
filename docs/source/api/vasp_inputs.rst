======================
vise.api.vasp_inputs
======================

.. module:: vise.api.vasp_inputs
   :synopsis: VASP input file generation

This module provides functions for generating VASP input files
(INCAR, KPOINTS, POSCAR, POTCAR) for various calculation types.

Quick Example
-------------

::

    from pymatgen.core import Structure
    from vise.api import vasp_inputs

    # Load structure
    struct = Structure.from_file("POSCAR")

    # Create VASP input set for structure optimization
    inputs = vasp_inputs.create_vasp_set(
        struct,
        task="structure_opt",  # or Task.structure_opt
        xc="pbe",              # or Xc.pbe
    )

    # Write input files to directory
    inputs.write("./relax")

    # Custom INCAR settings
    inputs = vasp_inputs.create_vasp_set(
        struct,
        task="band",
        xc="hse",
        user_incar_settings={"EDIFF": 1e-6, "NELM": 200}
    )
    inputs.write("./band")

Available Options
-----------------

Tasks
~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Task
     - Description
   * - ``structure_opt``
     - Structure optimization (relaxation)
   * - ``band``
     - Band structure calculation
   * - ``dos``
     - Density of states calculation
   * - ``dielectric_dfpt``
     - Dielectric function via DFPT
   * - ``dielectric_finite_field``
     - Dielectric function via finite field
   * - ``phonon_force``
     - Phonon force calculation
   * - ``defect``
     - Point defect calculation
   * - ``cluster_opt``
     - Molecule/cluster optimization

Exchange-Correlation Functionals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - XC
     - Description
   * - ``pbe``
     - PBE-GGA (standard DFT)
   * - ``pbesol``
     - PBEsol (solids optimized)
   * - ``lda``
     - Local density approximation
   * - ``scan``
     - SCAN meta-GGA
   * - ``hse``
     - HSE06 hybrid functional

Data Classes
------------

VaspInputSet
~~~~~~~~~~~~

.. autoclass:: VaspInputSet
   :members:
   :undoc-members:

Functions
---------

create_vasp_set
~~~~~~~~~~~~~~~

.. autofunction:: create_vasp_set

get_available_tasks
~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_available_tasks

get_available_xc
~~~~~~~~~~~~~~~~

.. autofunction:: get_available_xc

See Also
--------

- :doc:`structure` - Structure analysis before input generation
- :doc:`band` - Analyze band structure results
- CLI command: ``vise vi`` for VASP input generation
