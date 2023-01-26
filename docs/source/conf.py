# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import pathlib
import sys
import sphinx_rtd_theme
from pandemia.version import VERSION

sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())


project = 'pandemia'
copyright = '2023, Thompson, James AND Wattam, Stephen'
author = 'Thompson, James AND Wattam, Stephen'
release = VERSION

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx_rtd_theme',
    # 'myst_parser',
    'sphinx.ext.napoleon',
    'autoapi.extension',
    'm2r2',
]


templates_path = ['_templates']
exclude_patterns = []

autoapi_dirs = ['../../src/pandemia']
autoapi_type = "python"
autoapi_root = "api"

autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]

autoapi_keep_files = False
autodoc_typehints = "signature"
autoapi_add_toctree_entry = False

# source_suffix = {
#     '.rst': 'restructuredtext',
#     '.md': 'markdown',
# }

source_suffix = ['.rst', '.md']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"

html_static_path = ['_static']
