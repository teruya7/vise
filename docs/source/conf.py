# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'vise'
copyright = '2020, Yu Kumagai'
author = 'Yu Kumagai'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]

# Autodoc configuration
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'show-inheritance': True,
    'undoc-members': False,
}
autodoc_typehints = 'description'

# Napoleon settings (for Google/NumPy style docstrings)
napoleon_google_docstrings = True
napoleon_numpy_docstrings = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# Intersphinx mapping to external documentation
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'pymatgen': ('https://pymatgen.org/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# Theme configuration
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Theme options
html_theme_options = {
    # Navigation options
    'collapse_navigation': False,       # Don't collapse sidebar navigation
    'sticky_navigation': True,          # Keep navigation visible when scrolling
    'navigation_depth': 4,              # Show all levels in sidebar
    'includehidden': True,              # Include hidden toctrees
    'titles_only': False,               # Show all items, not just titles
    
    # Display options
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    
    # Style options
    'style_nav_header_background': '#0a1628',  # Deep sea darkest
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom CSS files
html_css_files = [
    'custom.css',
]

# Sidebar configuration - Show global TOC on all pages
html_sidebars = {
    '**': [
        'globaltoc.html',
        'searchbox.html',
    ]
}

# HTML context for templates
html_context = {
    'display_github': True,
    'github_user': 'kumagai-group',
    'github_repo': 'vise',
    'github_version': 'master',
    'conf_py_path': '/docs/source/',
}

# Logo (optional - uncomment if you have one)
# html_logo = '_static/logo.png'

# Favicon (optional - uncomment if you have one)
# html_favicon = '_static/favicon.ico'
