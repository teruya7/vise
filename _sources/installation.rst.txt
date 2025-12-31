Installation
============

This guide covers the installation of vise and its dependencies.

.. contents:: Table of Contents
   :local:
   :depth: 2

Requirements
------------

- Python 3.7 or higher
- VASP 5.4.4 or later (for running calculations)
- pymatgen POTCAR files (for input file generation)

Basic Installation
------------------

The latest stable version is available on PyPI:

.. code-block:: bash

    pip install vise

For the development version:

.. code-block:: bash

    git clone https://github.com/kumagai-group/vise.git
    cd vise
    pip install -e ".[dev]"

POTCAR Configuration
--------------------

vise uses pymatgen to generate POTCAR files. You must configure the path to your VASP pseudopotential files.

1. Create or edit ``~/.pmgrc.yaml``:

.. code-block:: yaml

    PMG_VASP_PSP_DIR: /path/to/your/potcar/directory

2. The directory structure should be:

::

    /path/to/your/potcar/directory/
    ├── POT_GGA_PAW_PBE/
    │   ├── Ac/
    │   │   └── POTCAR
    │   ├── Ag/
    │   │   └── POTCAR
    │   └── ...
    ├── POT_GGA_PAW_PBE_54/
    │   └── ...
    └── POT_LDA_PAW/
        └── ...

3. Verify the configuration:

.. code-block:: bash

    python -c "from pymatgen.io.vasp import Potcar; Potcar(['Si'])"

Optional Dependencies
---------------------

BoltzTraP2 (Effective Mass)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For effective mass calculations:

.. code-block:: bash

    pip install BoltzTraP2

Or install vise with the bolztrap extra:

.. code-block:: bash

    pip install vise[bolztrap]

Phonopy (Phonon Calculations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For phonon calculations:

.. code-block:: bash

    pip install phonopy

Materials Project API
~~~~~~~~~~~~~~~~~~~~~

To use Materials Project integration, add your API key to ``~/.pmgrc.yaml``:

.. code-block:: yaml

    PMG_MAPI_KEY: your_api_key_here

Get your API key at: https://materialsproject.org/api

Verification
------------

Verify your installation:

.. code-block:: bash

    # Check version
    vise --version

    # Check available commands
    vise --help

    # Test structure analysis (requires POSCAR)
    vise symmetry

Quick Test
~~~~~~~~~~

.. code-block:: python

    from pymatgen.core import Structure, Lattice
    from vise.api import vasp_inputs

    # Create a simple structure
    struct = Structure(
        Lattice.cubic(4.0),
        ["Si", "Si"],
        [[0, 0, 0], [0.25, 0.25, 0.25]]
    )

    # Generate VASP inputs
    inputs = vasp_inputs.create_vasp_set(struct, task="structure_opt", xc="pbe")
    print(inputs.incar)

Troubleshooting Installation
----------------------------

POTCAR Not Found
~~~~~~~~~~~~~~~~

If you see ``PmgVaspPspDirError``:

1. Check that ``PMG_VASP_PSP_DIR`` is set in ``~/.pmgrc.yaml``
2. Verify the directory path exists and contains POTCAR files
3. Check directory permissions

Import Errors
~~~~~~~~~~~~~

If you encounter import errors:

.. code-block:: bash

    # Reinstall with all dependencies
    pip install --force-reinstall vise

    # Or for development
    pip install -e ".[dev]" --force-reinstall

BoltzTraP2 Issues
~~~~~~~~~~~~~~~~~

BoltzTraP2 requires compiled extensions. On some systems:

.. code-block:: bash

    # Install build dependencies first
    pip install numpy scipy
    pip install BoltzTraP2
