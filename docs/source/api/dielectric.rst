===================
vise.api.dielectric
===================

.. module:: vise.api.dielectric
   :synopsis: Dielectric function analysis

This module provides functions for analyzing dielectric function
and optical properties from VASP output files.

Quick Example
-------------

::

    from vise.api import dielectric

    # Analyze dielectric function
    result = dielectric.analyze_dielectric("vasprun.xml", "OUTCAR")

    # Save plots for different components
    result.save_plot("absorption.pdf", plot_type="absorption_coeff")
    result.save_plot("real.pdf", plot_type="real")
    result.save_plot("imag.pdf", plot_type="imag")

    # Save data in multiple formats
    result.save_json("diele_func.json")
    result.save_csv("diele_func.csv")

    # Load from CSV for further analysis
    loaded = dielectric.load_dielectric_from_csv("diele_func.csv")

Plot Types
----------

The following plot types are available:

- ``absorption_coeff``: Absorption coefficient
- ``real``: Real part of dielectric function
- ``imag``: Imaginary part of dielectric function
- ``reflectivity``: Reflectivity
- ``refractive_index``: Refractive index
- ``extinction_coeff``: Extinction coefficient

Data Classes
------------

DielectricAnalysisResult
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: DielectricAnalysisResult
   :members:
   :undoc-members:

Functions
---------

analyze_dielectric
~~~~~~~~~~~~~~~~~~

.. autofunction:: analyze_dielectric

load_dielectric_from_csv
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: load_dielectric_from_csv

See Also
--------

- :doc:`band` - Band structure for electronic properties
- CLI command: ``vise pdfp`` for dielectric function plotting
