# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
"""VISE: VASP Integrated Supporting Environment.

VISE is a Python package that provides a high-level interface for
VASP (Vienna Ab-initio Simulation Package) calculations, including
input file generation, output analysis, and workflow automation.

This package exports the version string at the top level for easy access.
"""
from vise.version import __version__

__all__ = ["__version__"]
