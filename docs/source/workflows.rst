Workflow Examples
=================

This guide shows complete workflows for common calculation scenarios.

.. contents:: Table of Contents
   :local:
   :depth: 2

Structure Optimization Workflow
-------------------------------

A typical workflow for optimizing a crystal structure.

Step 1: Prepare Structure
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Option A: From Materials Project
    vise gp -f SrTiO3

    # Option B: From existing POSCAR
    # Just place your POSCAR file in the directory

Step 2: Check Symmetry
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    vise si

    # Output:
    # Space group: Pm-3m (221)
    # Crystal system: cubic
    # Point group: m-3m

Step 3: Generate Input Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Create input files
    vise vs -t structure_opt

    # Or with custom settings
    vise vs -t structure_opt -uis ENCUT 600 EDIFF 1e-7

Step 4: Run VASP
~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Run VASP (outside vise)
    mpirun -np 16 vasp_std

Step 5: Check Results
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Get band edge info
    vise be

    # Copy CONTCAR for next calculation
    cp CONTCAR ../band/POSCAR

Band Structure Calculation
--------------------------

Calculate and plot electronic band structure.

Step 1: Prepare from Relaxed Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    mkdir band && cd band
    cp ../relax/CONTCAR POSCAR

Step 2: Generate Band Calculation Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Get VBM/CBM from previous calculation
    vise be -v ../relax/vasprun.xml

    # Generate band inputs (automatically uses primitive cell)
    vise vs -t band --prev_dir ../relax

Step 3: Run VASP and Plot
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    mpirun -np 16 vasp_std

    # Plot band structure
    vise pb -y -5 5 -f band.pdf

    # With JSON data export
    vise pb --to-json

DOS Calculation
---------------

Calculate and plot density of states.

.. code-block:: bash

    mkdir dos && cd dos
    cp ../relax/CONTCAR POSCAR

    # Generate DOS inputs
    vise vs -t dos --prev_dir ../relax

    # Run VASP
    mpirun -np 16 vasp_std

    # Plot DOS
    vise pd -y -10 10 -f dos.pdf

HSE Calculation Workflow
------------------------

Hybrid functional calculations require more careful setup.

Step 1: PBE Pre-relaxation
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    mkdir pbe_relax && cd pbe_relax
    vise vs -t structure_opt
    mpirun -np 16 vasp_std

Step 2: HSE Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    mkdir ../hse_relax && cd ../hse_relax
    cp ../pbe_relax/CONTCAR POSCAR
    cp ../pbe_relax/WAVECAR .  # Speeds up convergence

    vise vs -t structure_opt -x hse
    mpirun -np 16 vasp_std

Step 3: HSE Band Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    mkdir ../hse_band && cd ../hse_band
    cp ../hse_relax/CONTCAR POSCAR
    cp ../hse_relax/WAVECAR .

    vise vs -t band -x hse --prev_dir ../hse_relax
    mpirun -np 16 vasp_std

    vise pb -f hse_band.pdf

Dielectric Function Calculation
-------------------------------

Calculate optical properties.

.. code-block:: bash

    mkdir diele && cd diele
    cp ../relax/CONTCAR POSCAR

    # Generate inputs for dielectric function
    vise vs -t dielectric_function --prev_dir ../relax

    mpirun -np 16 vasp_std

    # Plot absorption coefficient
    vise pdf -f absorption.pdf

    # Export to CSV
    vise pdf --to-csv

Phonon Calculation Workflow
---------------------------

Calculate phonon dispersion using phonopy.

Step 1: Create Supercell
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    mkdir phonon && cd phonon
    cp ../relax/CONTCAR POSCAR

    # Create 2x2x2 supercell
    vise mpp -s 2 2 2

Step 2: Run Force Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    vise vs -t phonon_force
    mpirun -np 16 vasp_std

Step 3: Create Phonon Band Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    vise mpf -f phonon_band.pdf

Python API Workflow
-------------------

For automation and high-throughput, use the Python API.

.. code-block:: python

    from pathlib import Path
    from pymatgen.core import Structure
    from vise.api import vasp_inputs, structure, band

    # Load and analyze structure
    struct = Structure.from_file("POSCAR")
    info = structure.get_symmetry_info(struct)
    print(f"Space group: {info.space_group_symbol}")

    # Get primitive cell for band calculation
    prim = structure.get_primitive(struct)

    # Generate inputs programmatically
    for task in ["structure_opt", "band", "dos"]:
        inputs = vasp_inputs.create_vasp_set(
            prim if task == "band" else struct,
            task=task,
            xc="pbe"
        )
        outdir = Path(task)
        inputs.write(outdir)
        print(f"Created inputs in {outdir}")

    # Analyze results
    result = band.analyze_band("band/vasprun.xml", "band/KPOINTS")
    print(f"Band gap: {result.band_gap:.3f} eV")

High-Throughput Example
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pathlib import Path
    from pymatgen.core import Structure
    from vise.api import vasp_inputs, materials_project

    # Get structures from Materials Project
    formulas = ["MgO", "CaO", "SrO", "BaO"]

    for formula in formulas:
        struct = materials_project.get_most_stable_structure(formula)

        # Create calculation directory
        calc_dir = Path(formula)
        calc_dir.mkdir(exist_ok=True)

        # Generate inputs
        inputs = vasp_inputs.create_vasp_set(
            struct,
            task="structure_opt",
            xc="pbe"
        )
        inputs.write(calc_dir)

        print(f"Created {formula} inputs")
