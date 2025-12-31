Tutorial -- vise.yaml
---------------------

Configuration files allow you to set default parameters for all vise operations.

.. contents:: Table of Contents
   :local:
   :depth: 2

Configuration File Location
===========================

vise searches for configuration files from the current directory up to the root:

- :code:`vise.yaml` (visible)
- :code:`.vise.yaml` (hidden)

Files in deeper directories override those in parent directories.

Example structure::

    ~/projects/
    ├── vise.yaml          # Project-wide defaults
    └── material_A/
        └── vise.yaml      # Material-specific overrides

Basic Configuration
===================

POTCAR Override
~~~~~~~~~~~~~~~

Override default POTCAR selections:

.. code-block:: yaml

    overridden_potcar:
      - Mg_pv
      - O_h

The Mg_pv and O_h POTCAR files will be used instead of default Mg and O.

XC Functional Default
~~~~~~~~~~~~~~~~~~~~~

Set default exchange-correlation functional:

.. code-block:: yaml

    xc: hse

INCAR Settings
~~~~~~~~~~~~~~

Override INCAR parameters:

.. code-block:: yaml

    user_incar_settings:
      LASPH: False
      ENCUT: 600
      EDIFF: 1e-8

Available Options
=================

All properties from the Defaults class can be set:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Key
     - Type
     - Description
   * - xc
     - string
     - Default XC functional (pbe, hse, scan, etc.)
   * - task
     - string
     - Default task type
   * - kpoint_density
     - float
     - K-point density in Å
   * - insulator_kpoint_density
     - float
     - K-point density for insulators
   * - defect_kpoint_density
     - float
     - K-point density for defect calculations
   * - symmetry_length_tolerance
     - float
     - Symmetry detection tolerance (Å)
   * - symmetry_angle_tolerance
     - float
     - Symmetry angle tolerance (degrees)
   * - str_opt_encut_factor
     - float
     - ENCUT multiplier for structure optimization
   * - overridden_potcar
     - list
     - Custom POTCAR selections
   * - user_incar_settings
     - dict
     - INCAR parameter overrides
   * - user_incar_tags
     - dict
     - Custom INCAR tag categories

Complete Example
================

.. code-block:: yaml

    # ~/projects/vise.yaml

    # Default functional for this project
    xc: pbesol

    # Higher k-point density
    kpoint_density: 3.0
    insulator_kpoint_density: 2.0

    # Custom POTCAR set
    overridden_potcar:
      - Zr_sv
      - O

    # Default INCAR overrides
    user_incar_settings:
      ENCUT: 600
      EDIFF: 1e-7
      ISMEAR: 0
      SIGMA: 0.05

    # Stricter symmetry tolerance
    symmetry_length_tolerance: 0.001

Priority Order
==============

Settings are applied in this order (later overrides earlier):

1. vise internal defaults
2. vise.yaml files from root to current directory
3. Command line arguments

Command line always has highest priority:

.. code-block:: bash

    # Overrides vise.yaml setting
    vise vs -x hse

See Also
========

- :doc:`installation` - POTCAR configuration
- :doc:`cli_reference` - Command line options
- :doc:`troubleshooting` - Configuration issues
