#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright (c) 2024. Distributed under the terms of the MIT License.

"""
Typer-based CLI for vise.

This module provides a clean, modern CLI using Typer.
It calls the vise.api module for all operations.
"""

import warnings
from pathlib import Path
from typing import List, Optional

import typer
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import UnknownPotcarWarning
from rich.console import Console
from rich.table import Table

from vise import __version__
from vise.defaults import defaults

# Suppress POTCAR warnings
warnings.simplefilter('ignore', UnknownPotcarWarning)

# Create main app
app = typer.Typer(
    name="vise",
    help="VASP Integrated Supporting Environment - Generate VASP inputs and analyze outputs.",
    add_completion=True,
    rich_markup_mode="rich",
)

console = Console()


def version_callback(value: bool):
    if value:
        console.print(f"vise version: {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None, "--version", "-V",
        callback=version_callback,
        is_eager=True,
        help="Show version and exit."
    ),
):
    """
    VISE - VASP Integrated Supporting Environment.
    
    Generate VASP input files and analyze calculation results.
    """
    pass


# =============================================================================
# structure_info command (alias: si)
# =============================================================================
@app.command("structure_info", rich_help_panel="Structure Analysis")
@app.command("si", hidden=True)
def structure_info(
    poscar: str = typer.Option(
        "POSCAR", "-p", "--poscar",
        help="POSCAR-type file name to read."
    ),
    symprec: float = typer.Option(
        defaults.symmetry_length_tolerance, "-s", "--symprec",
        help="Length tolerance in Å for symmetry analysis."
    ),
    angle_tolerance: float = typer.Option(
        defaults.symmetry_angle_tolerance, "-a", "--angle-tolerance",
        help="Angle tolerance in degrees for symmetry analysis."
    ),
    show_conventional: bool = typer.Option(
        False, "-c", "--conventional",
        help="Output the conventional cell."
    ),
    show_primitive: bool = typer.Option(
        False, "--primitive",
        help="Output the primitive cell."
    ),
):
    """
    Show structure symmetry information.
    
    Analyzes crystal structure and displays space group, crystal system,
    and lattice parameters.
    """
    from vise.api import structure

    struct = Structure.from_file(poscar)

    if show_primitive:
        primitive = structure.get_primitive(struct, symprec, angle_tolerance)
        if primitive != struct:
            typer.echo(primitive.to(fmt="poscar"))
        else:
            console.print("[yellow]Input structure is already primitive.[/yellow]")
    elif show_conventional:
        conventional = structure.get_conventional(struct, symprec, angle_tolerance)
        if conventional != struct:
            typer.echo(conventional.to(fmt="poscar"))
        else:
            console.print("[yellow]Input structure is already conventional.[/yellow]")
    else:
        info = structure.get_symmetry_info(struct, symprec, angle_tolerance)

        table = Table(title="Structure Information")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="green")

        table.add_row("Space group", f"{info.space_group_symbol} ({info.space_group_number})")
        table.add_row("Point group", info.point_group)
        table.add_row("Crystal system", info.crystal_system)
        table.add_row("Bravais lattice", info.bravais_lattice)
        table.add_row("Volume", f"{info.volume:.6f} Å³")
        table.add_row("a, b, c", f"{info.lattice_abc[0]:.6f}, {info.lattice_abc[1]:.6f}, {info.lattice_abc[2]:.6f} Å")
        table.add_row("α, β, γ", f"{info.lattice_angles[0]:.6f}, {info.lattice_angles[1]:.6f}, {info.lattice_angles[2]:.6f}°")
        table.add_row("Is primitive", "Yes" if info.is_primitive else "No")

        console.print(table)


# =============================================================================
# get_poscar command (alias: gp)
# =============================================================================
@app.command("get_poscar", rich_help_panel="Structure")
@app.command("gp", hidden=True)
def get_poscar(
    mpid: Optional[str] = typer.Option(
        None, "-m", "--mpid",
        help="Materials Project ID (e.g., mp-149)."
    ),
    formula: Optional[str] = typer.Option(
        None, "-f", "--formula",
        help="Chemical formula (e.g., SrTiO3). Gets the most stable structure."
    ),
):
    """
    Get POSCAR from Materials Project.
    
    Downloads structure from Materials Project database and saves as POSCAR.
    """
    from vise.api import materials_project

    if mpid is None and formula is None:
        console.print("[red]Error:[/red] Either --mpid or --formula must be specified.")
        raise typer.Exit(1)

    if mpid:
        console.print(f"Fetching structure for [cyan]{mpid}[/cyan]...")
        struct = materials_project.get_structure_by_id(mpid, save_poscar=True)
        info = materials_project.get_material_info(mpid)
        console.print("[green]✓[/green] Saved POSCAR")
        console.print(f"  Band gap: {info.band_gap:.3f} eV")
        console.print(f"  Magnetization: {info.total_magnetization:.3f} μB")
    else:
        console.print(f"Searching for [cyan]{formula}[/cyan]...")
        entries = materials_project.search_materials(formula)

        if not entries:
            console.print(f"[red]Error:[/red] No entries found for {formula}")
            raise typer.Exit(1)

        # Show candidates
        table = Table(title=f"Materials Project entries for {formula}")
        table.add_column("MP ID", style="cyan")
        table.add_column("E above hull", style="yellow")
        table.add_column("Space group")
        table.add_column("Band gap")

        for entry in entries[:10]:
            table.add_row(
                entry.material_id,
                f"{entry.energy_above_hull:.3f} eV/atom",
                entry.space_group,
                f"{entry.band_gap:.3f} eV"
            )
        console.print(table)

        # Get most stable
        struct = materials_project.get_structure_by_id(
            entries[0].material_id, save_poscar=True
        )
        console.print(f"[green]✓[/green] Saved POSCAR for {entries[0].material_id}")


# =============================================================================
# vasp_set command (alias: vs)
# =============================================================================
@app.command("vasp_set", rich_help_panel="VASP Input")
@app.command("vs", hidden=True)
def vasp_set(
    poscar: Optional[Path] = typer.Option(
        None, "-p", "--poscar",
        help="POSCAR file path. If not specified, searches in current directory."
    ),
    task: str = typer.Option(
        str(defaults.task), "-t", "--task",
        help="Task type: structure_opt, band, dos, etc."
    ),
    xc: str = typer.Option(
        str(defaults.xc), "-x", "--xc",
        help="Exchange-correlation functional: pbe, hse, scan, etc."
    ),
    kpt_density: Optional[float] = typer.Option(
        None, "-k", "--kpt-density",
        help="K-point density in Å. Default depends on task."
    ),
    potcar: Optional[List[str]] = typer.Option(
        None, "--potcar",
        help="Override POTCAR: e.g., Mg_pv O_h"
    ),
    user_incar_settings: Optional[List[str]] = typer.Option(
        None, "-uis", "--user-incar-settings",
        help="INCAR overrides: e.g., EDIFF 1e-6 NELM 200"
    ),
    output_dir: Path = typer.Option(
        Path.cwd(), "-d", "--dir",
        help="Output directory for input files."
    ),
):
    """
    Generate VASP input files.
    
    Creates INCAR, KPOINTS, POSCAR, and POTCAR for the specified task.
    """
    from vise.api import vasp_inputs

    # Load structure
    poscar_path = poscar or Path("POSCAR")
    if not poscar_path.exists():
        console.print(f"[red]Error:[/red] {poscar_path} not found.")
        raise typer.Exit(1)

    struct = Structure.from_file(str(poscar_path))

    # Parse user INCAR settings
    incar_settings = {}
    if user_incar_settings:
        it = iter(user_incar_settings)
        for key in it:
            try:
                value = next(it)
                # Try to convert to appropriate type
                try:
                    value = int(value)
                except ValueError:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                incar_settings[key] = value
            except StopIteration:
                console.print(f"[yellow]Warning:[/yellow] Odd number of INCAR settings, ignoring {key}")

    # Parse POTCAR overrides
    potcar_dict = None
    if potcar:
        potcar_dict = {}
        for p in potcar:
            if "_" in p:
                # e.g., "Mg_pv" -> {"Mg": "Mg_pv"}
                element = p.split("_")[0]
                potcar_dict[element] = p

    console.print("Generating VASP input files...")
    console.print(f"  Task: [cyan]{task}[/cyan]")
    console.print(f"  XC: [cyan]{xc}[/cyan]")

    try:
        inputs = vasp_inputs.create_vasp_set(
            struct,
            task=task,
            xc=xc,
            kpt_density=kpt_density,
            overridden_potcar=potcar_dict,
            user_incar_settings=incar_settings if incar_settings else None,
        )

        inputs.write(output_dir)
        console.print(f"[green]✓[/green] Created input files in {output_dir}")

        # Show INCAR summary
        table = Table(title="INCAR Summary")
        table.add_column("Tag", style="cyan")
        table.add_column("Value")

        important_tags = ["ENCUT", "EDIFF", "ISMEAR", "SIGMA", "IBRION", "NSW", "ALGO"]
        for tag in important_tags:
            if tag in inputs.incar:
                table.add_row(tag, str(inputs.incar[tag]))

        console.print(table)

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


# =============================================================================
# plot_band command (alias: pb)
# =============================================================================
@app.command("plot_band", rich_help_panel="Analysis")
@app.command("pb", hidden=True)
def plot_band(
    vasprun: Path = typer.Option(
        Path(defaults.vasprun), "-v", "--vasprun",
        help="Path to vasprun.xml file."
    ),
    kpoints: str = typer.Option(
        "KPOINTS", "-k", "--kpoints",
        help="Path to KPOINTS file."
    ),
    y_range: List[float] = typer.Option(
        [-10.0, 10.0], "-y", "--y-range",
        help="Energy range for plot."
    ),
    filename: str = typer.Option(
        "band.pdf", "-f", "--filename",
        help="Output PDF filename."
    ),
):
    """
    Plot band structure.
    
    Analyzes VASP band calculation results and creates a band structure plot.
    """
    from vise.api import band

    if not vasprun.exists():
        console.print(f"[red]Error:[/red] {vasprun} not found.")
        raise typer.Exit(1)

    console.print(f"Analyzing band structure from [cyan]{vasprun}[/cyan]...")

    result = band.analyze_band(str(vasprun), kpoints)

    # Show band edge info
    if result.is_metal:
        console.print("[yellow]System is metallic (no band gap)[/yellow]")
    else:
        console.print(f"Band gap: [green]{result.band_gap:.4f} eV[/green]")
        if result.vbm:
            console.print(f"VBM: {result.vbm:.4f} eV")
        if result.cbm:
            console.print(f"CBM: {result.cbm:.4f} eV")

    # Save plot
    result.save_plot(filename, energy_range=tuple(y_range))
    console.print(f"[green]✓[/green] Saved plot to {filename}")

    # Save JSON
    result.save_json()
    console.print("[green]✓[/green] Saved band_plot_info.json")


# =============================================================================
# plot_dos command (alias: pd)
# =============================================================================
@app.command("plot_dos", rich_help_panel="Analysis")
@app.command("pd", hidden=True)
def plot_dos(
    vasprun: Path = typer.Option(
        Path(defaults.vasprun), "-v", "--vasprun",
        help="Path to vasprun.xml file."
    ),
    outcar: Path = typer.Option(
        Path(defaults.outcar), "-o", "--outcar",
        help="Path to OUTCAR file."
    ),
    x_range: Optional[List[float]] = typer.Option(
        None, "-x", "--x-range",
        help="Energy range for plot."
    ),
    filename: str = typer.Option(
        "dos.pdf", "-f", "--filename",
        help="Output PDF filename."
    ),
    title: Optional[str] = typer.Option(
        None, "--title",
        help="Plot title."
    ),
):
    """
    Plot density of states.
    
    Analyzes VASP DOS calculation results and creates a DOS plot.
    """
    from vise.api import dos

    if not vasprun.exists():
        console.print(f"[red]Error:[/red] {vasprun} not found.")
        raise typer.Exit(1)

    console.print(f"Analyzing DOS from [cyan]{vasprun}[/cyan]...")

    result = dos.analyze_dos(str(vasprun), str(outcar))

    # Show band info
    if result.is_metal:
        console.print("[yellow]System is metallic[/yellow]")
        console.print(f"Fermi energy: {result.efermi:.4f} eV")
    else:
        console.print(f"Band gap: [green]{result.band_gap:.4f} eV[/green]")

    # Save plot
    energy_range = tuple(x_range) if x_range else (-10.0, 10.0)
    result.save_plot(filename, energy_range=energy_range, title=title)
    console.print(f"[green]✓[/green] Saved plot to {filename}")


# =============================================================================
# plot_diele_func command (alias: pdf)
# =============================================================================
@app.command("plot_diele_func", rich_help_panel="Analysis")
@app.command("pdf", hidden=True)
def plot_diele_func(
    vasprun: Path = typer.Option(
        Path(defaults.vasprun), "-v", "--vasprun",
        help="Path to vasprun.xml file."
    ),
    outcar: Path = typer.Option(
        Path(defaults.outcar), "-o", "--outcar",
        help="Path to OUTCAR file."
    ),
    plot_type: str = typer.Option(
        "absorption_coeff", "--plot-type",
        help="Plot type: real, imag, absorption_coeff"
    ),
    directions: List[str] = typer.Option(
        ["ave"], "-d", "--directions",
        help="Directions: xx, yy, zz, xy, yz, xz, ave"
    ),
    filename: Optional[str] = typer.Option(
        None, "-f", "--filename",
        help="Output PDF filename."
    ),
    calc_kk: bool = typer.Option(
        False, "--calc-kk",
        help="Calculate real part via Kramers-Kronig."
    ),
    to_csv: bool = typer.Option(
        False, "--to-csv",
        help="Export data to CSV."
    ),
    title: Optional[str] = typer.Option(
        None, "--title",
        help="Plot title."
    ),
):
    """
    Plot dielectric function.
    
    Plots dielectric function, absorption coefficient, or related quantities.
    """
    from vise.api import dielectric

    if not vasprun.exists():
        console.print(f"[red]Error:[/red] {vasprun} not found.")
        raise typer.Exit(1)

    console.print(f"Analyzing dielectric function from [cyan]{vasprun}[/cyan]...")

    result = dielectric.analyze_dielectric(
        str(vasprun), str(outcar),
        use_vasp_real=not calc_kk
    )

    # Save plot
    output_filename = filename or f"{plot_type}.pdf"
    result.save_plot(
        output_filename,
        plot_type=plot_type,
        directions=directions,
        title=title
    )
    console.print(f"[green]✓[/green] Saved plot to {output_filename}")

    # Save JSON
    result.save_json()
    console.print("[green]✓[/green] Saved diele_func_data.json")

    if to_csv:
        result.save_csv()
        console.print("[green]✓[/green] Saved CSV file")


# =============================================================================
# band_edge command (alias: be)
# =============================================================================
@app.command("band_edge", rich_help_panel="Analysis")
@app.command("be", hidden=True)
def band_edge(
    vasprun: Path = typer.Option(
        Path(defaults.vasprun), "-v", "--vasprun",
        help="Path to vasprun.xml file."
    ),
    outcar: Path = typer.Option(
        Path(defaults.outcar), "-o", "--outcar",
        help="Path to OUTCAR file."
    ),
):
    """
    Show band edge properties.
    
    Calculates and displays VBM, CBM, and band gap information.
    """
    from vise.api import band_edge as be_api

    if not vasprun.exists():
        console.print(f"[red]Error:[/red] {vasprun} not found.")
        raise typer.Exit(1)

    result = be_api.get_band_edge_properties(str(vasprun), str(outcar))

    if result.is_metal:
        console.print("[yellow]System is metallic (no band gap)[/yellow]")
        console.print(f"Fermi energy: {result.efermi:.4f} eV")
    else:
        table = Table(title="Band Edge Properties")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="green")

        gap_type = "direct" if result.is_direct else "indirect"
        table.add_row("Band gap", f"{result.band_gap:.4f} eV ({gap_type})")
        table.add_row("VBM energy", f"{result.vbm_info.energy:.4f} eV")
        table.add_row("VBM k-point", str(result.vbm_info.kpoint_coords))
        table.add_row("CBM energy", f"{result.cbm_info.energy:.4f} eV")
        table.add_row("CBM k-point", str(result.cbm_info.kpoint_coords))

        console.print(table)


# =============================================================================
# effective_mass command (alias: em)
# =============================================================================
@app.command("effective_mass", rich_help_panel="Analysis")
@app.command("em", hidden=True)
def effective_mass(
    vasprun: Path = typer.Option(
        Path(defaults.vasprun), "-v", "--vasprun",
        help="Path to vasprun.xml file."
    ),
    outcar: Path = typer.Option(
        Path(defaults.outcar), "-o", "--outcar",
        help="Path to OUTCAR file."
    ),
    temperature: float = typer.Option(
        300.0, "-t", "--temperature",
        help="Temperature in Kelvin."
    ),
    concentrations: Optional[List[float]] = typer.Option(
        None, "-c", "--concentrations",
        help="Carrier concentrations (as exponents, e.g., 18 19 20 for 1e18, 1e19, 1e20)."
    ),
):
    """
    Calculate effective mass.
    
    Uses BoltzTrap2 to calculate effective mass at given temperature and
    carrier concentrations.
    """
    from vise.api import band_edge as be_api

    if not vasprun.exists():
        console.print(f"[red]Error:[/red] {vasprun} not found.")
        raise typer.Exit(1)

    # Convert concentration exponents to actual values
    conc_list = None
    if concentrations:
        conc_list = [10 ** c for c in concentrations]

    console.print(f"Calculating effective mass at [cyan]{temperature} K[/cyan]...")

    try:
        result = be_api.calculate_effective_mass(
            str(vasprun), str(outcar),
            temperature=temperature,
            concentrations=conc_list
        )

        console.print(result)
        result.save_json()
        console.print("[green]✓[/green] Saved effective_mass.json")

    except ImportError:
        console.print("[red]Error:[/red] BoltzTrap2 is required for effective mass calculation.")
        console.print("Install with: pip install BoltzTrap2")
        raise typer.Exit(1)
    except ValueError as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


# =============================================================================
# UTILITY COMMANDS
# =============================================================================

# =============================================================================
# make_atom_poscars command (alias: map)
# =============================================================================
@app.command("make_atom_poscars", rich_help_panel="Utilities")
@app.command("map", hidden=True)
def make_atom_poscars(
    elements: Optional[List[str]] = typer.Option(
        None, "-e", "--elements",
        help="Element symbols (e.g., Si O Mg). If not specified, all elements.",
    ),
    output_dir: Path = typer.Option(
        Path.cwd(), "-d", "--dir",
        help="Directory where atom calculation directories are created.",
    ),
):
    """
    Create POSCAR files for isolated atom calculations.
    
    Generates directories with POSCARs for single-atom calculations,
    useful for reference energies.
    """
    from vise.api import util

    console.print("Creating atom POSCARs...")

    try:
        dirs = util.make_atom_poscars(elements=elements, output_dir=output_dir)
        console.print(f"[green]✓[/green] Created {len(dirs)} atom directories in {output_dir}")
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


# =============================================================================
# make_phonon_poscars command (alias: mpp)
# =============================================================================
@app.command("make_phonon_poscars", rich_help_panel="Utilities")
@app.command("mpp", hidden=True)
def make_phonon_poscars(
    unitcell: str = typer.Option(
        "POSCAR", "-u", "--unitcell",
        help="Path to unitcell POSCAR file.",
    ),
    supercell_matrix: List[int] = typer.Option(
        ..., "-s", "--supercell",
        help="Supercell matrix (3 integers for diagonal, e.g., 2 2 2).",
    ),
    output_dir: Path = typer.Option(
        Path.cwd(), "-d", "--dir",
        help="Output directory.",
    ),
):
    """
    Create phonon calculation setup.
    
    Generates a supercell POSCAR and phonopy_input.json for phonon calculations.
    """
    from vise.api import util

    if not Path(unitcell).exists():
        console.print(f"[red]Error:[/red] {unitcell} not found.")
        raise typer.Exit(1)

    console.print(f"Creating phonon setup with supercell matrix {supercell_matrix}...")

    try:
        result = util.create_phonon_setup(
            unitcell=unitcell,
            supercell_matrix=supercell_matrix,
            output_dir=output_dir
        )
        console.print(f"[green]✓[/green] Created POSCAR ({len(result.supercell)} atoms)")
        console.print(f"[green]✓[/green] Created {result.phonopy_input_file}")
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


# =============================================================================
# make_phonon_figs command (alias: mpf)
# =============================================================================
@app.command("make_phonon_figs", rich_help_panel="Utilities")
@app.command("mpf", hidden=True)
def make_phonon_figs(
    phonopy_input: str = typer.Option(
        "phonopy_input.json", "-pi", "--phonopy-input",
        help="Path to phonopy_input.json file.",
    ),
    vasprun: str = typer.Option(
        "vasprun.xml", "-vn", "--vasprun",
        help="Path to vasprun.xml from force calculation.",
    ),
    filename: str = typer.Option(
        "phonon_band.pdf", "-f", "--filename",
        help="Output plot filename.",
    ),
):
    """
    Create phonon band structure plot.
    
    Analyzes phonon calculation results and creates a band structure plot.
    """
    from vise.api import util

    if not Path(phonopy_input).exists():
        console.print(f"[red]Error:[/red] {phonopy_input} not found.")
        raise typer.Exit(1)

    if not Path(vasprun).exists():
        console.print(f"[red]Error:[/red] {vasprun} not found.")
        raise typer.Exit(1)

    console.print("Analyzing phonon calculation...")

    try:
        result = util.analyze_phonon(phonopy_input, vasprun, filename)
        console.print(f"[green]✓[/green] Saved plot to {result}")
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


# =============================================================================
# spin_decomposed_volumetric_files command (alias: sdvf)
# =============================================================================
@app.command("spin_decomposed_volumetric_files", rich_help_panel="Utilities")
@app.command("sdvf", hidden=True)
def spin_decomposed_volumetric_files(
    chgcar: str = typer.Option(
        ..., "-c", "--chgcar",
        help="CHGCAR-type file name with spin data.",
    ),
):
    """
    Create spin-decomposed volumetric files.
    
    Splits a spin-polarized CHGCAR into spin-up and spin-down components.
    """
    from vise.api import util

    if not Path(chgcar).exists():
        console.print(f"[red]Error:[/red] {chgcar} not found.")
        raise typer.Exit(1)

    console.print(f"Creating spin-decomposed files from [cyan]{chgcar}[/cyan]...")

    try:
        up_file, down_file = util.make_spin_decomposed_volumetric_files(chgcar)
        console.print(f"[green]✓[/green] Created {up_file}")
        console.print(f"[green]✓[/green] Created {down_file}")
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


# =============================================================================
# light_weight_vol_data command (alias: lwvd)
# =============================================================================
@app.command("light_weight_vol_data", rich_help_panel="Utilities")
@app.command("lwvd", hidden=True)
def light_weight_vol_data(
    volumetric_file: str = typer.Option(
        ..., "-v", "--volumetric-file",
        help="CHGCAR-type volumetric file name.",
    ),
    output_filename: Optional[str] = typer.Option(
        None, "-o", "--output",
        help="Created file name. Default: {input}_lw",
    ),
    output_vesta: Optional[str] = typer.Option(
        None, "--output-vesta",
        help="Output VESTA file name.",
    ),
    original_vesta: Optional[str] = typer.Option(
        None, "--original-vesta",
        help="Original VESTA file to use as base.",
    ),
):
    """
    Create lightweight volumetric data file.
    
    Reduces file size while preserving isosurface information.
    Optionally creates a VESTA file with isosurfaces.
    """
    from vise.api import util

    if not Path(volumetric_file).exists():
        console.print(f"[red]Error:[/red] {volumetric_file} not found.")
        raise typer.Exit(1)

    console.print(f"Creating lightweight volumetric data from [cyan]{volumetric_file}[/cyan]...")

    try:
        lw_file = util.make_light_weight_volumetric_data(
            volumetric_file,
            output_filename=output_filename
        )
        console.print(f"[green]✓[/green] Created {lw_file}")

        if output_vesta:
            vesta_file = util.create_vesta_file(
                volumetric_file,
                output_vesta,
                original_vesta_file=original_vesta,
                lw_volumetric_file=lw_file
            )
            console.print(f"[green]✓[/green] Created {vesta_file}")
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


# =============================================================================
# Entry point
# =============================================================================
if __name__ == "__main__":
    app()

