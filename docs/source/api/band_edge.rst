==================
vise.api.band_edge
==================

.. module:: vise.api.band_edge
   :synopsis: Band edge properties analysis

This module provides functions for analyzing band edge properties
including VBM/CBM details and effective mass calculations.

Quick Example
-------------

Band Edge Properties
~~~~~~~~~~~~~~~~~~~~

::

    from vise.api import band_edge

    # Get band edge properties
    result = band_edge.get_band_edge_properties("vasprun.xml", "OUTCAR")

    if result.is_metal:
        print(f"Metallic (Fermi: {result.efermi} eV)")
    else:
        print(f"Band gap: {result.band_gap:.4f} eV")
        print(f"Type: {'direct' if result.is_direct else 'indirect'}")
        print(f"VBM: {result.vbm_info.energy:.4f} eV at {result.vbm_info.kpoint_coords}")
        print(f"CBM: {result.cbm_info.energy:.4f} eV at {result.cbm_info.kpoint_coords}")

Effective Mass Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Requires BoltzTrap2 to be installed.

::

    from vise.api import band_edge

    # Calculate effective mass
    result = band_edge.calculate_effective_mass(
        "vasprun.xml", "OUTCAR",
        temperature=300.0,
        concentrations=[1e18, 1e19, 1e20]
    )

    print(result)
    result.save_json("effective_mass.json")

Data Classes
------------

BandEdgeResult
~~~~~~~~~~~~~~

.. autoclass:: BandEdgeResult
   :members:
   :undoc-members:

EffectiveMassResult
~~~~~~~~~~~~~~~~~~~

.. autoclass:: EffectiveMassResult
   :members:
   :undoc-members:

Functions
---------

get_band_edge_properties
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_band_edge_properties

calculate_effective_mass
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: calculate_effective_mass

Dependencies
------------

.. note::

   Effective mass calculation requires `BoltzTrap2 <https://gitlab.com/sousaw/BoltzTraP2>`_
   to be installed::

       pip install BoltzTraP2

See Also
--------

- :doc:`band` - Full band structure analysis
- :doc:`dos` - Density of states analysis
