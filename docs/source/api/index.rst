===============
API Reference
===============

This section provides comprehensive API documentation for vise, automatically
generated from source code docstrings.

.. contents:: Modules
   :local:
   :depth: 1

Overview
--------

The ``vise.api`` module provides a clean, high-level Python interface for:

- **Structure Analysis**: Symmetry detection, primitive/conventional cell conversion
- **VASP Input Generation**: Automated INCAR, KPOINTS, POSCAR, POTCAR creation
- **Band Structure Analysis**: Band gap, VBM/CBM analysis, and plotting
- **DOS Analysis**: Density of states extraction and visualization
- **Dielectric Function**: Optical property analysis
- **Band Edge Properties**: Effective mass, carrier concentration
- **Materials Project Integration**: Structure search and download

Quick Import
------------

All API functions are accessible via the ``vise.api`` submodules::

    from vise.api import structure, vasp_inputs, band, dos
    from vise.api import dielectric, band_edge, materials_project, util

Module Reference
----------------

.. toctree::
   :maxdepth: 2
   
   structure
   vasp_inputs
   band
   dos
   dielectric
   band_edge
   materials_project
   util

Legacy API (Low-level)
----------------------

For more fine-grained control, you can use the lower-level API directly.

::

    from pathlib import Path
    
    from pymatgen.core import Structure, Lattice
    
    from vise.input_set.input_options import CategorizedInputOptions
    from vise.input_set.vasp_input_files import VaspInputFiles
    from vise.input_set.task import Task
    from vise.input_set.xc import Xc
    
    structure = Structure(Lattice.cubic(1), ["Mg", "O"], [[0.0]*3, [0.5]*3])
    
    categorized_input_options = CategorizedInputOptions(
                structure=structure,
                task=Task.band,
                xc=Xc.pbe, 
                overridden_potcar={"Mg": "Mg_pv"})
    
    input_files = VaspInputFiles(categorized_input_options, overridden_incar_settings={"NSW": 20})
    input_files.create_input_files(dirname=Path("."))

The CategorizedInputOptions class constructor takes the keyword arguments that 
are arguments of 
`generate_potcar <https://github.com/kumagai-group/vise/blob/master/vise/input_set/potcar_generator.py>`_,
`IncarSettingsGenerator <https://github.com/kumagai-group/vise/blob/master/vise/input_set/incar_settings_generator.py>`_,
and `StructureKpointsGenerator <https://github.com/kumagai-group/vise/blob/master/vise/input_set/structure_kpoints_generator.py>`_.
An example is overridden_potcar shown above.

The VaspInputFiles class constructor also takes the overridden_incar_settings, which can control the INCAR tags.

Note also that the vise.yaml files are also parsed.

FireWorks Integration
---------------------

Here, we show an example of FireTask in `FireWorks <https://materialsproject.github.io/fireworks/index.html>`_.


::

    from pathlib import Path
    
    from fireworks import FiretaskBase, explicit_serialize
    from vise.input_set.input_options import CategorizedInputOptions
    from vise.input_set.vasp_input_files import VaspInputFiles
    
    @explicit_serialize
    class WriteVaspInputsTask(FiretaskBase):
    
        required_params = ["task", "xc"]
        optional_params = ["input_options", "overridden_incar_settings"]
    
        def run_task(self, fw_spec):
            categorized_input_options = CategorizedInputOptions(
                structure=fw_spec["structure"],
                task=self["task"],
                xc=self["xc"],
                **self.get("input_options", {}))
    
            input_files = VaspInputFiles(categorized_input_options,
                                         self.get("overridden_incar_settings", {}))
            input_files.create_input_files(dirname=Path("."))
