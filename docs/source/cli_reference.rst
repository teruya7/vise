CLI Reference
=============

Complete reference for all vise command-line commands.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
--------

vise provides a modern CLI built with Typer. All commands have short aliases for convenience.

.. code-block:: bash

    # Show version
    vise --version

    # Show help
    vise --help

    # Get help for a specific command
    vise vasp_set --help

Command Summary
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Command
     - Alias
     - Description
   * - structure_info
     - si
     - Show structure symmetry information
   * - get_poscar
     - gp
     - Get POSCAR from Materials Project
   * - vasp_set
     - vs
     - Generate VASP input files
   * - plot_band
     - pb
     - Plot band structure
   * - plot_dos
     - pd
     - Plot density of states
   * - plot_diele_func
     - pdf
     - Plot dielectric function
   * - band_edge
     - be
     - Show band edge properties
   * - effective_mass
     - em
     - Calculate effective mass
   * - make_atom_poscars
     - map
     - Create atom POSCARs for reference
   * - make_phonon_poscars
     - mpp
     - Create phonon calculation setup
   * - make_phonon_figs
     - mpf
     - Create phonon band plots
   * - spin_decomposed_volumetric_files
     - sdvf
     - Create spin-decomposed volumetric files
   * - light_weight_vol_data
     - lwvd
     - Create lightweight volumetric data

Input Generation Commands
-------------------------

vasp_set (vs)
~~~~~~~~~~~~~

Generate VASP input files (INCAR, KPOINTS, POSCAR, POTCAR).

**Options:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Short
     - Description
   * - ``--poscar``
     - ``-p``
     - Path to POSCAR file
   * - ``--task``
     - ``-t``
     - Calculation task (structure_opt, band, dos, etc.)
   * - ``--xc``
     - ``-x``
     - XC functional (pbe, pbesol, hse, scan, etc.)
   * - ``--user-incar-settings``
     - ``-uis``
     - INCAR overrides (key value pairs)
   * - ``--dir``
     - ``-d``
     - Output directory

**Examples:**

.. code-block:: bash

    # Basic structure optimization
    vise vs

    # HSE band structure calculation
    vise vs -t band -x hse

    # With custom INCAR settings
    vise vs -t dos -uis EDIFF 1e-6 NELM 200

    # Output to specific directory
    vise vs -t structure_opt -d ./relax

**Available Tasks:**

- ``structure_opt`` - Full structure optimization
- ``structure_opt_tight`` - High-precision optimization
- ``structure_opt_rough`` - Quick rough optimization
- ``cluster_opt`` - Cluster/molecule optimization
- ``band`` - Band structure calculation
- ``dos`` - Density of states
- ``defect`` - Defect calculation
- ``defect_2d`` - 2D defect calculation
- ``dielectric_dfpt`` - Dielectric via DFPT
- ``dielectric_finite_field`` - Dielectric via finite field
- ``dielectric_function`` - Optical dielectric function
- ``phonon_force`` - Phonon force calculation

**Available XC Functionals:**

- ``pbe`` - PBE (default)
- ``pbesol`` - PBEsol
- ``lda`` - LDA
- ``scan`` - SCAN meta-GGA
- ``hse`` - HSE06 hybrid
- ``pbe0`` - PBE0 hybrid

Structure Commands
------------------

structure_info (si)
~~~~~~~~~~~~~~~~~~~

Show structure symmetry information.

**Options:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Short
     - Description
   * - ``--poscar``
     - ``-p``
     - Path to POSCAR file (default: POSCAR)
   * - ``--symprec``
     - ``-s``
     - Symmetry precision in Å (default: 0.01)
   * - ``--conventional``
     - ``-c``
     - Output conventional cell
   * - ``--primitive``
     -
     - Output primitive cell

**Examples:**

.. code-block:: bash

    # Show symmetry info
    vise si

    # With custom tolerance
    vise si -s 0.1

    # Get primitive cell
    vise si --primitive

get_poscar (gp)
~~~~~~~~~~~~~~~

Download structure from Materials Project.

**Options:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Short
     - Description
   * - ``--mpid``
     - ``-m``
     - Materials Project ID (e.g., mp-149)
   * - ``--formula``
     - ``-f``
     - Chemical formula (gets most stable)

**Examples:**

.. code-block:: bash

    # Get by Materials Project ID
    vise gp -m mp-149

    # Get most stable structure by formula
    vise gp -f SrTiO3

Analysis Commands
-----------------

plot_band (pb)
~~~~~~~~~~~~~~

Plot band structure from VASP output.

**Options:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Short
     - Description
   * - ``--vasprun``
     - ``-v``
     - Path to vasprun.xml
   * - ``--kpoints``
     - ``-k``
     - Path to KPOINTS file
   * - ``--y-range``
     - ``-y``
     - Energy range (eV)
   * - ``--filename``
     - ``-f``
     - Output filename (default: band.pdf)

**Examples:**

.. code-block:: bash

    # Basic plot
    vise pb

    # Custom energy range
    vise pb -y -5 5

    # Custom output file
    vise pb -f my_band.png

plot_dos (pd)
~~~~~~~~~~~~~

Plot density of states.

**Options:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Short
     - Description
   * - ``--vasprun``
     - ``-v``
     - Path to vasprun.xml
   * - ``--outcar``
     - ``-o``
     - Path to OUTCAR
   * - ``--y-range``
     - ``-y``
     - Energy range (eV)
   * - ``--filename``
     - ``-f``
     - Output filename (default: dos.pdf)
   * - ``--title``
     -
     - Plot title

**Examples:**

.. code-block:: bash

    vise pd
    vise pd -y -10 10 -f dos.png

plot_diele_func (pdf)
~~~~~~~~~~~~~~~~~~~~~

Plot dielectric function and related properties.

**Options:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Short
     - Description
   * - ``--vasprun``
     - ``-v``
     - Path to vasprun.xml
   * - ``--outcar``
     - ``-o``
     - Path to OUTCAR
   * - ``--y-range``
     - ``-y``
     - Energy range (eV)
   * - ``--calc-kk``
     - ``-ckk``
     - Use Kramers-Kronig transformation
   * - ``--to-csv``
     -
     - Export to CSV
   * - ``--filename``
     - ``-f``
     - Output filename

**Examples:**

.. code-block:: bash

    # Plot absorption coefficient
    vise pdf

    # Use Kramers-Kronig transformation
    vise pdf -ckk

    # Export data to CSV
    vise pdf --to-csv

band_edge (be)
~~~~~~~~~~~~~~

Show band edge properties (VBM, CBM, band gap).

**Examples:**

.. code-block:: bash

    vise be
    vise be -v path/to/vasprun.xml

effective_mass (em)
~~~~~~~~~~~~~~~~~~~

Calculate effective mass using BoltzTraP2.

**Options:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Short
     - Description
   * - ``--vasprun``
     - ``-v``
     - Path to vasprun.xml
   * - ``--outcar``
     - ``-o``
     - Path to OUTCAR
   * - ``--temperature``
     - ``-temp``
     - Temperature in K (default: 300)
   * - ``--concentrations``
     - ``-c``
     - Carrier concentrations (exponents)

**Examples:**

.. code-block:: bash

    # Default temperature, typical concentrations
    vise em -c 17 18 19 20

    # At 500K
    vise em -temp 500 -c 18 19

Utility Commands
----------------

make_atom_poscars (map)
~~~~~~~~~~~~~~~~~~~~~~~

Create POSCARs for isolated atom reference calculations.

**Examples:**

.. code-block:: bash

    # All elements
    vise map

    # Specific elements
    vise map --elements Si O Mg

make_phonon_poscars (mpp)
~~~~~~~~~~~~~~~~~~~~~~~~~

Create phonon calculation setup with supercell.

**Examples:**

.. code-block:: bash

    # 2x2x2 supercell
    vise mpp -s 2 2 2

    # From specific unitcell
    vise mpp -u POSCAR_prim -s 3 3 3

make_phonon_figs (mpf)
~~~~~~~~~~~~~~~~~~~~~~

Create phonon band structure plot from calculated forces.

**Examples:**

.. code-block:: bash

    vise mpf

spin_decomposed_volumetric_files (sdvf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Split spin-polarized volumetric data into components.

**Examples:**

.. code-block:: bash

    vise sdvf -c CHGCAR

light_weight_vol_data (lwvd)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create lightweight volumetric data for visualization.

**Examples:**

.. code-block:: bash

    vise lwvd -c CHGCAR
    vise lwvd -c CHGCAR --output-vesta result.vesta
