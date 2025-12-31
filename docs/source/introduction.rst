=====================
Introduction of vise
=====================

First-principles calculations are becoming increasingly prevalent in materials science.
However, preparing input files for the calculations and analyzing calculation
results are still not an easy task for researchers other than experts,
and human errors could come up.

To overcome the situation, we have developed :code:`vise` code,
which means *VASP Integrated Supporting Environment*.
:code:`Vise` supports :code:`VASP` users to generate its input files
for several tasks with suitable defaults, and allows the users to
tune some parameters depending on their purposes.

Key Features
------------

**Input Generation**

- Automatic INCAR parameter selection based on calculation task
- Smart k-point mesh generation proportional to reciprocal lattice
- POTCAR selection with customizable presets

**Analysis Tools**

- Electronic band structure plotting
- Density of states (DOS) analysis
- Dielectric function and optical absorption
- Effective mass calculation (via BoltzTraP2)
- Band edge property extraction

**Integration**

- Materials Project database access
- Phonon calculation support (via phonopy)
- Simple Python API for automation
- Modern CLI with Typer

Usage Modes
-----------

**Command Line Interface**

.. code-block:: bash

    # Generate VASP inputs
    vise vs -t structure_opt

    # Plot band structure
    vise pb

    # Get structure from Materials Project
    vise gp -m mp-149

See :doc:`cli_reference` for complete command documentation.

**Python API**

.. code-block:: python

    from vise.api import vasp_inputs, band

    # Generate inputs programmatically
    inputs = vasp_inputs.create_vasp_set(structure, task="band")
    inputs.write("./band")

    # Analyze results
    result = band.analyze_band("vasprun.xml", "KPOINTS")

See :doc:`vise_api` for complete API documentation.

Quick Start
-----------

1. Install vise: ``pip install vise``
2. Configure POTCAR: See :doc:`installation`
3. Follow the tutorials:
   - :doc:`tutorial_preparation_poscar`
   - :doc:`tutorial_input_set`
   - :doc:`tutorial_properties`

**Note: Units used in vise are eV for energy and Å for length
following the VASP convention.**

Next Steps
----------

- :doc:`installation` - Detailed installation instructions
- :doc:`workflows` - Step-by-step workflow examples
- :doc:`troubleshooting` - Common issues and solutions
