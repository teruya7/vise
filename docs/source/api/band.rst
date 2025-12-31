==============
vise.api.band
==============

.. module:: vise.api.band
   :synopsis: Band structure analysis

This module provides functions for analyzing band structure
from VASP output files.

Quick Example
-------------

::

    from vise.api import band

    # Analyze band structure from VASP output
    result = band.analyze_band("vasprun.xml", "KPOINTS")

    # Check if metallic
    if result.is_metal:
        print("System is metallic")
    else:
        print(f"Band gap: {result.band_gap:.4f} eV")
        print(f"VBM: {result.vbm:.4f} eV")
        print(f"CBM: {result.cbm:.4f} eV")

    # Save plot
    result.save_plot("band.pdf", energy_range=(-5, 5))

    # Save data as JSON
    result.save_json("band_data.json")

Data Classes
------------

BandAnalysisResult
~~~~~~~~~~~~~~~~~~

.. autoclass:: BandAnalysisResult
   :members:
   :undoc-members:

Functions
---------

analyze_band
~~~~~~~~~~~~

.. autofunction:: analyze_band

Output Files
------------

The analysis can generate the following output files:

- **band.pdf** (or .png): Publication-quality band structure plot
- **band_plot_info.json**: Complete band data in JSON format for further processing

See Also
--------

- :doc:`dos` - Density of states analysis
- :doc:`band_edge` - Band edge properties (VBM/CBM details)
- CLI command: ``vise pb`` for band plotting
