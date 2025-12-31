==========================
vise.api.materials_project
==========================

.. module:: vise.api.materials_project
   :synopsis: Materials Project integration

This module provides functions for searching and retrieving structures
from the Materials Project database.

Prerequisites
-------------

.. important::

   This module requires a Materials Project API key configured in ``~/.pmgrc.yaml``::

       PMG_MAPI_KEY: your_api_key_here

   Get your API key at: https://next-gen.materialsproject.org/api

Quick Example
-------------

::

    from vise.api import materials_project

    # Search for materials by formula
    entries = materials_project.search_materials("TiO2")
    for entry in entries[:5]:
        print(f"{entry.material_id}: {entry.space_group}, Eg={entry.band_gap:.2f} eV")

    # Get structure by Materials Project ID
    struct = materials_project.get_structure_by_id("mp-149", save_poscar=True)

    # Get most stable structure for a formula
    struct = materials_project.get_most_stable_structure("LiFePO4", save_poscar=True)

    # Get detailed material info
    info = materials_project.get_material_info("mp-149")
    print(f"Band gap: {info.band_gap} eV")
    print(f"Magnetization: {info.total_magnetization} μB")

Data Classes
------------

MaterialEntry
~~~~~~~~~~~~~

.. autoclass:: MaterialEntry
   :members:
   :undoc-members:

MaterialInfo
~~~~~~~~~~~~

.. autoclass:: MaterialInfo
   :members:
   :undoc-members:

Functions
---------

search_materials
~~~~~~~~~~~~~~~~

.. autofunction:: search_materials

get_structure_by_id
~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_structure_by_id

get_most_stable_structure
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_most_stable_structure

get_material_info
~~~~~~~~~~~~~~~~~

.. autofunction:: get_material_info

See Also
--------

- :doc:`structure` - Analyze retrieved structures
- :doc:`vasp_inputs` - Generate VASP inputs from MP structures
- `Materials Project <https://materialsproject.org>`_ - Official website
