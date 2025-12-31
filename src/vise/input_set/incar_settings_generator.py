# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""INCAR settings generation for VASP calculations.

This module provides classes for generating INCAR settings based on the
calculation task type, exchange-correlation functional, and various
user-specified options.

Key Components:
    - IncarSettingsGenerator: Main class that generates complete INCAR settings
    - TaskIncarSettings: Generates task-specific INCAR parameters
    - XcIncarSettings: Generates XC functional-specific INCAR parameters
"""
from __future__ import annotations

from math import ceil
from typing import Dict, List, Optional, Union

from pymatgen.core import Composition, Structure
from pymatgen.io.vasp.sets import Potcar

from vise.analyzer.band_edge_properties import is_band_gap
from vise.defaults import defaults
from vise.input_set.datasets.dataset_util import LDAU, calc_kpar, num_bands
from vise.input_set.fft_grids import vasp_grid
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger

logger = get_logger(__name__)

# Type alias for INCAR setting values
IncarValue = Union[str, bool, int, float, List[int], List[float]]


class IncarSettingsGenerator:
    """Generator for complete INCAR settings based on calculation parameters.

    This class generates appropriate INCAR parameters based on the structure,
    task type, exchange-correlation functional, and various computational
    options. It handles settings for:
    - Default VASP parameters (LASPH, NELM, SIGMA, etc.)
    - Task-specific settings (PREC, LREAL, EDIFF, ISIF, etc.)
    - XC functional settings (ALGO, GGA, METAGGA, etc.)
    - Optional features (spin-orbit, Hubbard U, spectrum calculations)

    Attributes:
        incar_settings: Dictionary containing the generated INCAR parameters.
    """

    def __init__(
        self,
        structure: Structure,
        symbol_list: List[str],
        num_kpts: int,
        num_kpt_multiplication_factor: int,
        potcar: Potcar,
        task: Task,
        xc: Xc,
        dos_step_size: float = defaults.dos_step_size,
        charge: float = 0.0,
        band_gap: Optional[float] = None,
        vbm_cbm: Optional[List[float]] = None,
        exchange_ratio: float = 0.25,
        set_hubbard_u: Optional[bool] = None,
        auto_kpar: bool = True,
        cutoff_energy: Optional[float] = None,
        is_magnetization: bool = False,
        num_cores_for_kpar: int = defaults.default_num_cores,
        unused_core_ratio_threshold: float = defaults.unused_core_ratio_threshold,
        str_opt_encut_multi_factor: float = defaults.str_opt_encut_factor,
        multiples_for_grids: Optional[List[int]] = None,
        set_spin_orbit: bool = False,
    ) -> None:
        """Initialize the INCAR settings generator.

        Args:
            structure: pymatgen Structure object.
            symbol_list: List of element symbols in the structure.
            num_kpts: Number of k-points in the mesh.
            num_kpt_multiplication_factor: K-point multiplication factor for
                hybrid functional calculations.
            potcar: pymatgen Potcar object.
            task: Calculation task type (e.g., structure_opt, band, dos).
            xc: Exchange-correlation functional (e.g., pbe, hse).
            dos_step_size: Energy step size for DOS calculations in eV.
            charge: Net charge of the system (positive = cation).
            band_gap: Band gap in eV (None for metals).
            vbm_cbm: List of [VBM, CBM] in absolute eV.
            exchange_ratio: Exact exchange ratio for hybrid functionals.
            set_hubbard_u: Whether to apply Hubbard U corrections.
            auto_kpar: Whether to automatically set KPAR.
            cutoff_energy: Manual ENCUT override in eV.
            is_magnetization: Whether to perform spin-polarized calculation.
            num_cores_for_kpar: Number of CPU cores for KPAR calculation.
            unused_core_ratio_threshold: Threshold for unused core ratio in KPAR.
            str_opt_encut_multi_factor: ENCUT multiplier for structure optimization.
            multiples_for_grids: FFT grid multiples [NGX, NGY, NGZ].
            set_spin_orbit: Whether to enable spin-orbit coupling.
        """
        self._composition: Composition = structure.composition
        self._lattice = structure.lattice
        self._symbol_list = symbol_list
        self._num_kpt_multiplication_factor = num_kpt_multiplication_factor
        self._potcar = potcar
        self._num_kpts = num_kpts
        self._task = task
        self._xc = xc
        self._dos_step_size = dos_step_size
        self._charge = charge
        self._band_gap = band_gap
        self._vbm_cbm = vbm_cbm
        self._exchange_ratio = exchange_ratio
        self._auto_kpar = auto_kpar
        self._cutoff_energy = cutoff_energy
        self._is_magnetization = is_magnetization
        self._num_cores_for_kpar = num_cores_for_kpar
        self._unused_core_ratio_threshold = unused_core_ratio_threshold
        self._str_opt_encut_multi_factor = str_opt_encut_multi_factor
        self._multiples_for_grids = multiples_for_grids
        self._set_spin_orbit = set_spin_orbit

        self._incar_settings: Dict[str, IncarValue] = {}
        self._set_incar_settings(set_hubbard_u)

    def _set_incar_settings(self, set_hubbard_u: Optional[bool]) -> None:
        """Build complete INCAR settings based on all configuration options.

        Args:
            set_hubbard_u: Whether to apply Hubbard U corrections.
        """
        self._set_default_settings()
        self._set_task_related_settings()
        self._set_xc_related_settings()
        self._set_options_related_settings()
        self._set_spin_orbit_related_settings()

        if self._task.is_spectrum_task:
            self._set_spectrum_related_settings()
        if self._task.is_dielectric:
            self._set_dielectric_related_settings()
        if self._xc.is_hybrid_functional:
            self._set_hybrid_func_related_settings()
        if self._need_hubbard_u(set_hubbard_u):
            self._set_hubbard_u_related_settings()
        if self._auto_kpar:
            self._set_kpar()

        self._remove_incar_settings_with_none_values()

    def _set_default_settings(self) -> None:
        """Set default INCAR parameters common to all calculations."""
        self._incar_settings.update({
            "LASPH": True,    # Use aspherical PAW corrections
            "NELM": 100,      # Maximum electronic iterations
            "SIGMA": 0.1,     # Smearing width in eV
            "LCHARG": False,  # Don't write CHGCAR by default
        })

    def _set_task_related_settings(self) -> None:
        """Set INCAR parameters based on calculation task type."""
        task_settings = TaskIncarSettings(self._task)
        self._incar_settings.update({
            "PREC": task_settings.prec,
            "LREAL": task_settings.lreal,
            "EDIFF": task_settings.ediff,
            "ADDGRID": task_settings.addgrid_optional,
            "ISIF": task_settings.isif,
            "IBRION": task_settings.ibrion,
            "EDIFFG": task_settings.ediffg_optional,
            "NSW": task_settings.nsw,
            "POTIM": task_settings.potim_optional,
            "LORBIT": task_settings.lorbit,
        })

    def _set_xc_related_settings(self) -> None:
        """Set INCAR parameters based on exchange-correlation functional."""
        xc_settings = XcIncarSettings(self._xc)
        self._incar_settings.update({
            "ALGO": xc_settings.algo,
            "LWAVE": xc_settings.lwave,
            "GGA": xc_settings.gga_optional,
            "METAGGA": xc_settings.metagga_optional,
        })

    def _set_options_related_settings(self) -> None:
        """Set INCAR parameters based on user-provided options."""
        self._incar_settings.update({
            "ENCUT": self._encut,
            "ISMEAR": self._ismear,
            "ISPIN": self._ispin,
            "NBANDS": self._nbands,
            "NELECT": self._nelect,
        })

        if self._multiples_for_grids:
            self._set_fft_grid_settings()

    def _set_fft_grid_settings(self) -> None:
        """Set custom FFT grid dimensions based on multiples."""
        grids = [
            vasp_grid(
                self._incar_settings["ENCUT"],
                length,
                self._incar_settings["PREC"],
            )
            for length in self._lattice.abc
        ]
        multiples = self._multiples_for_grids
        ngx = ceil(grids[0] / multiples[0]) * multiples[0]
        ngy = ceil(grids[1] / multiples[1]) * multiples[1]
        ngz = ceil(grids[2] / multiples[2]) * multiples[2]
        self._incar_settings.update({"NGX": ngx, "NGY": ngy, "NGZ": ngz})

    def _set_spin_orbit_related_settings(self) -> None:
        """Enable spin-orbit coupling if requested."""
        if self._set_spin_orbit:
            self._incar_settings["LSORBIT"] = True

    def _need_hubbard_u(self, set_hubbard_u: Optional[bool]) -> bool:
        """Determine if Hubbard U corrections should be applied.

        Args:
            set_hubbard_u: User-specified Hubbard U setting.

        Returns:
            True if Hubbard U corrections should be applied.
        """
        if isinstance(set_hubbard_u, bool):
            return set_hubbard_u
        return False

    @property
    def _encut(self) -> float:
        """Calculate the plane-wave cutoff energy.

        Returns:
            ENCUT value in eV, either user-specified or derived from POTCAR.
        """
        if self._cutoff_energy:
            return self._cutoff_energy

        max_enmax = max([p.enmax for p in self._potcar])
        if self._task.is_lattice_relaxed:
            max_enmax *= self._str_opt_encut_multi_factor

        return round(max_enmax, 3)

    @property
    def _ismear(self) -> int:
        """Determine the smearing method.

        Returns:
            ISMEAR value: 0 (Gaussian), -4/-5 (tetrahedron methods).
        """
        # Tetrahedron method fails for irreducible k-points NKPT < 4 in VASP
        if self._task is Task.band:
            return 0

        if is_band_gap(self._band_gap, self._vbm_cbm, show_info=True) and self._num_kpts >= 4:
            if self._task in (Task.dos, Task.dielectric_function):
                # -4 and -5 show same results for spectra in tests
                return -4
            else:
                return -5

        return 0

    @property
    def _ispin(self) -> Optional[int]:
        """Determine spin polarization setting.

        Returns:
            ISPIN value: 2 for spin-polarized, None for non-spin-polarized.
        """
        if self._is_magnetization or self._task.need_spin:
            return 2
        return None

    @property
    def _nbands(self) -> Optional[int]:
        """Calculate number of bands for spectrum calculations.

        Returns:
            NBANDS value for plot tasks, None otherwise.
        """
        if self._task.is_plot_task:
            return num_bands(self._composition, self._potcar, self._set_spin_orbit)
        return None

    @property
    def _nelect(self) -> Optional[float]:
        """Calculate total electron count for charged systems.

        Returns:
            NELECT value if charge is non-zero, None otherwise.
        """
        if self._charge:
            element_composition = self._composition.element_composition
            nelect = sum([
                element_composition[potcar_single.element] * potcar_single.ZVAL
                for potcar_single in self._potcar
            ])
            return nelect - self._charge
        return None

    def _set_spectrum_related_settings(self) -> None:
        """Set INCAR parameters for DOS/spectrum calculations."""
        # Need one more step for VASP to remove weird huge values
        if self._vbm_cbm:
            emin = ceil(self._vbm_cbm[0]) - 15 - self._dos_step_size
            emax = ceil(self._vbm_cbm[1]) + 15
        else:
            emin = -20 - self._dos_step_size
            emax = 20

        self._incar_settings["EMIN"] = emin
        self._incar_settings["EMAX"] = emax
        num_steps = round((emax - emin) / self._dos_step_size)
        self._incar_settings["NEDOS"] = num_steps + 1

    def _set_dielectric_related_settings(self) -> None:
        """Set INCAR parameters for dielectric property calculations."""
        if self._task == Task.dielectric_dfpt:
            self._incar_settings["LEPSILON"] = True
        elif self._task == Task.dielectric_finite_field:
            self._incar_settings["LCALCEPS"] = True
            self._incar_settings["POTIM"] = 0.015
        elif self._task == Task.dielectric_function:
            self._incar_settings["LOPTICS"] = True
            # Real part of dielectric function depends on VASP version,
            # so we use CSHIFT=0.0 and apply KK transformation using vise
            self._incar_settings["CSHIFT"] = 0.0

    def _set_hybrid_func_related_settings(self) -> None:
        """Set INCAR parameters for hybrid functional calculations."""
        self._incar_settings["LHFCALC"] = True
        self._incar_settings["PRECFOCK"] = "Fast"
        self._incar_settings["TIME"] = 0.4
        self._incar_settings["AEXX"] = self._exchange_ratio

        if self._num_kpt_multiplication_factor != 1:
            self._incar_settings["NKRED"] = self._num_kpt_multiplication_factor

        if self._xc is Xc.hse:
            self._incar_settings["HFSCREEN"] = 0.208

    def _set_hubbard_u_related_settings(self) -> None:
        """Set INCAR parameters for DFT+U calculations."""
        ldau = LDAU(self._symbol_list)
        if ldau.is_ldau_needed:
            self._incar_settings["LDAU"] = True
            self._incar_settings["LDAUTYPE"] = 2
            self._incar_settings["LMAXMIX"] = ldau.lmaxmix
            self._incar_settings["LDAUPRINT"] = 1
            self._incar_settings["LDAUU"] = ldau.ldauu
            self._incar_settings["LDAUL"] = ldau.ldaul

    def _set_kpar(self) -> None:
        """Set KPAR for k-point parallelization."""
        if self._kpar_incompatible():
            return

        self._incar_settings["KPAR"] = calc_kpar(
            self._num_kpts,
            self._num_cores_for_kpar,
            self._unused_core_ratio_threshold,
        )

    def _kpar_incompatible(self) -> bool:
        """Check if KPAR is incompatible with other settings.

        Returns:
            True if KPAR cannot be used (e.g., with LELF=True).
        """
        if self._incar_settings.get("LELF", False):
            logger.warning(
                "Because KPAR cannot be set with LELF=True, KPAR is switched off."
            )
            return True
        return False

    def _remove_incar_settings_with_none_values(self) -> None:
        """Remove INCAR settings with None values."""
        for tag_name in list(self._incar_settings.keys()):
            if self._incar_settings[tag_name] is None:
                self._incar_settings.pop(tag_name)

    @property
    def incar_settings(self) -> Dict[str, IncarValue]:
        """Return the generated INCAR settings dictionary.

        Returns:
            Dictionary mapping INCAR tag names to their values.
        """
        return self._incar_settings


class TaskIncarSettings:
    """Generator for task-specific INCAR parameters.

    This class provides INCAR parameter values appropriate for different
    calculation types (structure optimization, band structure, DOS, etc.).

    Attributes:
        isif: Stress tensor calculation mode.
        ediff: Electronic convergence criterion.
        ediffg_optional: Ionic convergence criterion (optional).
        ibrion: Ionic relaxation algorithm.
        lorbit: Orbital projection method.
        lreal: Real-space projection setting.
        prec: Calculation precision.
        nsw: Maximum ionic steps.
        potim_optional: Ionic step size (optional).
        addgrid_optional: Additional FFT grid flag (optional).
    """

    def __init__(self, task: Task) -> None:
        """Initialize with the calculation task type.

        Args:
            task: The calculation task type.
        """
        self._task = task

    @property
    def isif(self) -> int:
        """Determine stress tensor calculation mode.

        Returns:
            ISIF value: 3 (full relaxation), 2 (ion only), 0 (no relaxation).

        Raises:
            NotImplementedError: If task type is not recognized.
        """
        if self._task.is_lattice_relaxed:
            return 3
        elif self._task.is_atom_relaxed_lattice_fixed or self._task is Task.phonon_force:
            return 2
        elif self._task in (Task.band, Task.dos) or self._task.is_dielectric:
            return 0
        else:
            raise NotImplementedError

    @property
    def ediff(self) -> float:
        """Determine electronic convergence criterion.

        Note: During dielectric_dfpt calculations, EDIFF is tightened automatically.

        Returns:
            EDIFF value in eV.

        Raises:
            NotImplementedError: If task type is not recognized.
        """
        if self._task.is_tight_calc:
            return 1e-8
        elif self._task in (Task.structure_opt, Task.cluster_opt):
            return 1e-7
        elif self._task in (Task.dielectric_dfpt, Task.dielectric_finite_field):
            return 1e-6
        elif self._task in (
            Task.dielectric_function,
            Task.band,
            Task.dos,
            Task.defect,
            Task.defect_2d,
        ):
            return 1e-5
        elif self._task in (Task.structure_opt_rough,):
            return 1e-4
        else:
            raise NotImplementedError

    @property
    def ediffg_optional(self) -> Optional[float]:
        """Determine ionic convergence criterion.

        Returns:
            EDIFFG value (negative = force-based) or None if not applicable.

        Raises:
            NotImplementedError: If task type is not recognized for relaxation.
        """
        if self._task.is_atom_relaxed:
            if self._task is Task.structure_opt_tight:
                return -0.001
            elif self._task in (Task.structure_opt, Task.cluster_opt):
                return -0.005
            elif self._task in (Task.defect, Task.defect_2d):
                return -0.03
            elif self._task is Task.structure_opt_rough:
                return -0.2
            else:
                raise NotImplementedError
        return None

    @property
    def ibrion(self) -> int:
        """Determine ionic relaxation algorithm.

        Returns:
            IBRION value: 8 (DFPT) or 2 (conjugate gradient).
        """
        return 8 if self._task is Task.dielectric_dfpt else 2

    @property
    def lorbit(self) -> Optional[int]:
        """Determine orbital projection method.

        Returns:
            LORBIT value: 10, 11, or None for phonon_force.
        """
        if self._task is Task.phonon_force:
            return None
        return 10 if self._task is not Task.dos else 11

    @property
    def lreal(self) -> Union[str, bool]:
        """Determine real-space projection setting.

        Returns:
            'Auto' for defect calculations, False otherwise.
        """
        return "Auto" if self._task in (Task.defect, Task.defect_2d) else False

    @property
    def prec(self) -> str:
        """Determine calculation precision.

        Returns:
            'Accurate' for tight calculations, 'Normal' otherwise.
        """
        return "Accurate" if self._task.is_tight_calc else "Normal"

    @property
    def nsw(self) -> int:
        """Determine maximum ionic relaxation steps.

        Returns:
            20 for relaxation tasks, 1 for single-point calculations.
        """
        return 20 if self._task.is_atom_relaxed else 1

    @property
    def potim_optional(self) -> Optional[float]:
        """Determine ionic step size for rough optimization.

        Returns:
            POTIM value or None if default should be used.
        """
        if self._task is Task.structure_opt_rough:
            return 0.1
        return None

    @property
    def addgrid_optional(self) -> Optional[bool]:
        """Determine if additional FFT grid is needed.

        Returns:
            True for tight calculations, None otherwise.
        """
        if self._task.is_tight_calc:
            return True
        return None


class XcIncarSettings:
    """Generator for exchange-correlation functional INCAR parameters.

    This class provides INCAR parameter values appropriate for different
    XC functionals (LDA, GGA, meta-GGA, hybrids).

    Attributes:
        algo: Electronic minimization algorithm.
        lwave: Whether to write WAVECAR.
        gga_optional: GGA type (optional).
        metagga_optional: Meta-GGA type (optional).
    """

    def __init__(self, xc: Xc) -> None:
        """Initialize with the XC functional type.

        Args:
            xc: The exchange-correlation functional.
        """
        self._xc = xc

    @property
    def algo(self) -> str:
        """Determine electronic minimization algorithm.

        Returns:
            'Damped' for hybrid functionals, 'Normal' otherwise.
        """
        return "Damped" if self._xc.is_hybrid_functional else "Normal"

    @property
    def lwave(self) -> bool:
        """Determine if WAVECAR should be written.

        Returns:
            True for hybrid functionals (needed for restart), False otherwise.
        """
        return True if self._xc.is_hybrid_functional else False

    @property
    def gga_optional(self) -> Optional[str]:
        """Determine GGA type for non-standard functionals.

        Returns:
            'PS' for PBEsol, None otherwise.
        """
        return "PS" if self._xc is Xc.pbesol else None

    @property
    def metagga_optional(self) -> Optional[str]:
        """Determine meta-GGA type.

        Returns:
            Meta-GGA name in uppercase, or None if not meta-GGA.
        """
        if self._xc.is_metagga:
            return self._xc.value.upper()
        return None
