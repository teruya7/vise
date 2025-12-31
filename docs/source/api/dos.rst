=============
vise.api.dos
=============

.. module:: vise.api.dos
   :synopsis: Density of states analysis

This module provides functions for analyzing density of states (DOS)
from VASP output files.

Quick Example
-------------

::

    from vise.api import dos

    # Analyze DOS from VASP output
    result = dos.analyze_dos("vasprun.xml", "OUTCAR")

    # Get basic information
    print(f"Band gap: {result.band_gap} eV")
    print(f"Fermi energy: {result.efermi} eV")

    # Save plot
    result.save_plot("dos.pdf", energy_range=(-10, 10))

    # Save data
    result.save_json("dos_data.json")

Data Classes
------------

DosAnalysisResult
~~~~~~~~~~~~~~~~~

.. autoclass:: DosAnalysisResult
   :members:
   :undoc-members:

Functions
---------

analyze_dos
~~~~~~~~~~~

.. autofunction:: analyze_dos

Output Files
------------

The analysis can generate the following output files:

- **dos.pdf** (or .png): Publication-quality DOS plot
- **dos_data.json**: Complete DOS data in JSON format

See Also
--------

- :doc:`band` - Band structure analysis
- :doc:`band_edge` - Band edge properties
- CLI command: ``vise pd`` for DOS plotting
