Change Log
==========

Version 0.9.5
-------------

**Release Date:** December 2024

Improvements
~~~~~~~~~~~~

- **Code Quality:** Comprehensive docstrings and type hints added throughout the codebase
- **Documentation:** Added installation guide, CLI reference, workflow examples, and troubleshooting guide

Changes
~~~~~~~

- Refactored input_set module for improved readability
- Refactored analyzer module with better type annotations
- Enhanced defaults.py with comprehensive property documentation

Version 0.9.4
-------------

API Modernization
~~~~~~~~~~~~~~~~~

- Added high-level ``vise.api`` module for easy Python access
- Migrated CLI from argparse to Typer for improved user experience
- Modernized project structure to src-layout

Version 0.9.0
-------------

Initial Features
~~~~~~~~~~~~~~~~

- VASP input file generation (INCAR, KPOINTS, POSCAR, POTCAR)
- Band structure analysis and plotting
- DOS analysis and plotting
- Dielectric function analysis
- Materials Project integration
- vise.yaml configuration system

Supported Tasks
~~~~~~~~~~~~~~~

- Structure optimization (standard, tight, rough)
- Band structure calculations
- DOS calculations
- Dielectric property calculations (DFPT, finite field, optical)
- Phonon force calculations
- Defect calculations

Supported XC Functionals
~~~~~~~~~~~~~~~~~~~~~~~~

- PBE, PBEsol, LDA
- SCAN meta-GGA
- HSE06, PBE0 hybrid functionals
