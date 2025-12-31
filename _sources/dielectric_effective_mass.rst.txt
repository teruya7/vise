Dielectric Function & Effective Mass
=====================================

Analyze optical properties and calculate effective mass from VASP calculations.

.. contents:: Table of Contents
   :local:
   :depth: 2

Dielectric Function
-------------------

Command Line
~~~~~~~~~~~~

.. code-block:: bash

    # Plot absorption coefficient
    vise pdf

    # Use Kramers-Kronig transformation
    vise pdf -ckk

    # Export to CSV
    vise pdf --to-csv

The ``-ckk`` option calculates the real part of the dielectric function 
from the imaginary part using Kramers-Kronig transformation.

Python API
~~~~~~~~~~

.. code-block:: python

    from vise.api import dielectric

    result = dielectric.analyze_dielectric("vasprun.xml", "OUTCAR")
    result.save_plot("absorption.pdf", plot_type="absorption_coeff")
    result.save_csv("diele_func.csv")

Effective Mass
--------------

Requires BoltzTraP2 installation. See :doc:`installation`.

Command Line
~~~~~~~~~~~~

.. code-block:: bash

    vise em -c 17 18 19 20

The ``-c`` option specifies carrier concentrations as exponents of base 10 in cm⁻³.

Python API
~~~~~~~~~~

.. code-block:: python

    from vise.api import band_edge

    result = band_edge.calculate_effective_mass(
        "vasprun.xml", "OUTCAR",
        temperature=300.0,
        concentrations=[1e17, 1e18, 1e19, 1e20]
    )
    result.save_json("effective_mass.json")

See Also
--------

- :doc:`band_dos` - Band structure and DOS analysis
- :doc:`cli_commands` - Full CLI reference
