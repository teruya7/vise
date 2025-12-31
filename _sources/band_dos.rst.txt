Band Structure & DOS
====================

Analyze electronic band structure and density of states from VASP calculations.

.. contents:: Table of Contents
   :local:
   :depth: 2

Band Structure
--------------

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

Density of States
-----------------

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

Band Edge Properties
--------------------

.. code-block:: bash

    vise be

.. code-block:: python

    from vise.api import band_edge

    result = band_edge.get_band_edge_properties("vasprun.xml", "OUTCAR")

    if not result.is_metal:
        print(f"Band gap: {result.band_gap:.4f} eV")
        print(f"Type: {'direct' if result.is_direct else 'indirect'}")

See Also
--------

- :doc:`dielectric_effective_mass` - Optical properties and effective mass
- :doc:`cli_commands` - Full CLI reference
