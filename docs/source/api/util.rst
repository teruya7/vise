==============
vise.api.util
==============

.. module:: vise.api.util
   :synopsis: Utility functions

This module provides various utility functions for common computational
materials science tasks.

Quick Example
-------------

Atom POSCARs
~~~~~~~~~~~~

Create isolated atom POSCARs for reference energy calculations::

    from vise.api import util

    # Create atom POSCARs for specific elements
    util.make_atom_poscars(["Si", "O"], output_dir="./atoms")

Phonon Calculations
~~~~~~~~~~~~~~~~~~~

::

    from vise.api import util

    # Create phonon calculation setup
    setup = util.create_phonon_setup("POSCAR", [2, 2, 2], output_dir="./phonon")
    print(f"Supercell atoms: {len(setup.supercell)}")

    # After VASP calculation, analyze phonon results
    util.analyze_phonon("phonopy_input.json", "vasprun.xml", "phonon_band.pdf")

Volumetric Data
~~~~~~~~~~~~~~~

::

    from vise.api import util

    # Create spin-decomposed volumetric files
    up, down = util.make_spin_decomposed_volumetric_files("CHGCAR")

    # Create lightweight volumetric data (for large files)
    util.make_light_weight_volumetric_data("CHGCAR")

Functions
---------

Atom Calculations
~~~~~~~~~~~~~~~~~

.. autofunction:: make_atom_poscars

Phonon Calculations
~~~~~~~~~~~~~~~~~~~

.. autofunction:: create_phonon_setup

.. autofunction:: analyze_phonon

Volumetric Data Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: make_spin_decomposed_volumetric_files

.. autofunction:: make_light_weight_volumetric_data

Data Classes
------------

PhonopySetup
~~~~~~~~~~~~

.. autoclass:: PhonopySetup
   :members:
   :undoc-members:

Dependencies
------------

.. note::

   Phonon-related functions require `phonopy <https://phonopy.github.io/phonopy/>`_
   to be installed::

       pip install phonopy

See Also
--------

- :doc:`vasp_inputs` - Generate VASP inputs for calculations
- :doc:`structure` - Structure analysis utilities
