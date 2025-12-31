Troubleshooting
===============

Common issues and solutions when using vise.

.. contents:: Table of Contents
   :local:
   :depth: 2

Installation Issues
-------------------

POTCAR Not Found
~~~~~~~~~~~~~~~~

**Error:**

.. code-block::

    PmgVaspPspDirError: Set PMG_VASP_PSP_DIR=<directory-path> in .pmgrc.yaml

**Solution:**

1. Create ``~/.pmgrc.yaml`` with:

   .. code-block:: yaml

       PMG_VASP_PSP_DIR: /path/to/your/potcar/directory

2. Verify directory structure:

   ::

       /path/to/your/potcar/directory/
       └── POT_GGA_PAW_PBE_54/
           ├── Ac/POTCAR
           ├── Ag/POTCAR
           └── ...

3. Test:

   .. code-block:: bash

       python -c "from pymatgen.io.vasp import Potcar; print(Potcar(['Si']))"

BoltzTraP2 Import Error
~~~~~~~~~~~~~~~~~~~~~~~

**Error:**

.. code-block::

    ImportError: BoltzTraP2 has to be installed and working

**Solution:**

1. Install BoltzTraP2:

   .. code-block:: bash

       pip install BoltzTraP2

2. If compilation fails, install build dependencies first:

   .. code-block:: bash

       pip install numpy scipy cython
       pip install BoltzTraP2

Materials Project API Key
~~~~~~~~~~~~~~~~~~~~~~~~~

**Error:**

.. code-block::

    MPRestError: API key not found

**Solution:**

Add your API key to ``~/.pmgrc.yaml``:

.. code-block:: yaml

    PMG_MAPI_KEY: your_api_key_here

Get your key at: https://materialsproject.org/api

Input Generation Issues
-----------------------

Wrong POTCAR Selected
~~~~~~~~~~~~~~~~~~~~~

**Problem:** vise selects unexpected POTCAR variant.

**Solution:** Override POTCAR selection:

.. code-block:: bash

    vise vs --potcar Mg_pv O

Or in Python:

.. code-block:: python

    inputs = vasp_inputs.create_vasp_set(
        struct,
        overridden_potcar={"Mg": "Mg_pv"}
    )

KPOINTS Not Gamma-Centered
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem:** Need Gamma-centered k-points for specific calculations.

**Solution:** Use ``--options``:

.. code-block:: bash

    vise vs --options gamma_centered True

Symmetry Detection Issues
~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem:** Wrong space group detected.

**Solution:** Adjust symmetry tolerance:

.. code-block:: bash

    vise si -s 0.1  # Increase tolerance

    # Or in vise.yaml
    symmetry_length_tolerance: 0.1
    symmetry_angle_tolerance: 10.0

Analysis Issues
---------------

Band Gap Shows Zero
~~~~~~~~~~~~~~~~~~~

**Problem:** Non-zero band gap system shows band_gap=0.

**Possible causes:**

1. **Metallic k-point sampling:** Need denser k-mesh for insulators

   .. code-block:: bash

       vise vs -t band --options kpt_density 2.0

2. **Spin-polarized calculation needed:** Add magnetization

   .. code-block:: bash

       vise vs -uis ISPIN 2

Band Structure Plot Missing Symmetry Points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem:** Band plot doesn't show high-symmetry labels.

**Solution:** Ensure KPOINTS file is in line mode:

.. code-block:: bash

    # Use band task, not dos
    vise vs -t band

DOS Shows Unexpected Gaps
~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem:** DOS has strange gaps or spikes.

**Possible causes:**

1. Too few k-points - increase density
2. Smearing too small - use SIGMA=0.1 or higher
3. Need more NEDOS points

VASP Calculation Issues
-----------------------

VASP Crashes with KPAR
~~~~~~~~~~~~~~~~~~~~~~

**Problem:** VASP crashes when KPAR is set.

**Solution:** Disable automatic KPAR:

.. code-block:: bash

    vise vs -uis KPAR 1

Or in vise.yaml:

.. code-block:: yaml

    options:
      auto_kpar: False

Slow HSE Calculations
~~~~~~~~~~~~~~~~~~~~~

**Tips for faster HSE:**

1. Start from PBE WAVECAR:

   .. code-block:: bash

       cp pbe_calc/WAVECAR hse_calc/

2. Use appropriate NKRED for k-point reduction:

   .. code-block:: bash

       vise vs -x hse --options num_kpt_factor 2

3. Consider using PRECFOCK = Fast (default in vise)

Memory Issues
~~~~~~~~~~~~~

**Problem:** Out of memory during calculation.

**Solutions:**

1. Reduce NBANDS:

   .. code-block:: bash

       vise vs -uis NBANDS 100

2. Use LREAL = Auto for large cells:

   .. code-block:: bash

       vise vs -t defect  # Automatic for defect task

3. Reduce k-point density:

   .. code-block:: bash

       vise vs --options kpt_density 3.0

Configuration Issues
--------------------

vise.yaml Not Being Read
~~~~~~~~~~~~~~~~~~~~~~~~

**Problem:** Settings in vise.yaml are ignored.

**Checklist:**

1. File name is exactly ``vise.yaml`` or ``.vise.yaml``
2. File is in current directory or parent directories
3. YAML syntax is valid:

   .. code-block:: bash

       python -c "import yaml; yaml.safe_load(open('vise.yaml'))"

Priority Order
~~~~~~~~~~~~~~

Settings are applied in this order (later overrides earlier):

1. vise defaults
2. vise.yaml in parent directories (root → current)
3. vise.yaml in current directory
4. Command line arguments

Getting Help
------------

1. Check command help:

   .. code-block:: bash

       vise <command> --help

2. View API documentation:

   .. code-block:: python

       from vise.api import vasp_inputs
       help(vasp_inputs.create_vasp_set)

3. Report issues: https://github.com/kumagai-group/vise/issues
