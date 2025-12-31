Tutorial -- Calculations of Various Properties
----------------------------------------------

This tutorial covers analyzing VASP calculation results.

.. contents:: Table of Contents
   :local:
   :depth: 2

Supported Properties
====================

.. csv-table:: Properties
   :file: properties.csv
   :header-rows: 1

**Notes:**

- Band paths are determined using `seekpath <https://www.materialscloud.org/work/tools/seekpath>`_. Please cite the related paper if used in publications.
- For carrier concentration, DOS or absorption coefficient calculations are required.

Band Structure Analysis
=======================

Command Line
~~~~~~~~~~~~

.. code-block:: bash

    # Basic plot
    vise pb

    # Custom energy range
    vise pb -y -5 5

    # Custom output file
    vise pb -f band.png

Python API
~~~~~~~~~~

.. code-block:: python

    from vise.api import band

    result = band.analyze_band("vasprun.xml", "KPOINTS")

    if result.is_metal:
        print("Metallic")
    else:
        print(f"Band gap: {result.band_gap:.4f} eV")
        print(f"VBM: {result.vbm:.4f} eV")
        print(f"CBM: {result.cbm:.4f} eV")

    result.save_plot("band.pdf", energy_range=(-5, 5))
    result.save_json("band_data.json")

DOS Analysis
============

Command Line
~~~~~~~~~~~~

.. code-block:: bash

    vise pd
    vise pd -y -10 10 -f dos.png

Python API
~~~~~~~~~~

.. code-block:: python

    from vise.api import dos

    result = dos.analyze_dos("vasprun.xml", "OUTCAR")
    print(f"Band gap: {result.band_gap} eV")
    print(f"Fermi energy: {result.efermi} eV")

    result.save_plot("dos.pdf", energy_range=(-10, 10))

Dielectric Function and Absorption
==================================

Command Line
~~~~~~~~~~~~

.. code-block:: bash

    # Plot absorption coefficient
    vise pdf

    # Use Kramers-Kronig transformation
    vise pdf -ckk

    # Export to CSV
    vise pdf --to-csv

The ``-ckk`` option calculates the real part of the dielectric function from the imaginary part using Kramers-Kronig transformation. The complex shift η is set via ``--ita``.

Python API
~~~~~~~~~~

.. code-block:: python

    from vise.api import dielectric

    result = dielectric.analyze_dielectric("vasprun.xml", "OUTCAR")
    result.save_plot("absorption.pdf", plot_type="absorption_coeff")
    result.save_csv("diele_func.csv")

Band Edge Properties
====================

Command Line
~~~~~~~~~~~~

.. code-block:: bash

    vise be

Python API
~~~~~~~~~~

.. code-block:: python

    from vise.api import band_edge

    result = band_edge.get_band_edge_properties("vasprun.xml", "OUTCAR")

    if result.is_metal:
        print(f"Metallic (Fermi: {result.efermi} eV)")
    else:
        print(f"Band gap: {result.band_gap:.4f} eV")
        print(f"Type: {'direct' if result.is_direct else 'indirect'}")

Effective Mass
==============

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
========

- :doc:`workflows` - Complete calculation workflows
- :doc:`cli_reference` - Full CLI command reference
- :doc:`vise_api` - Complete Python API reference
