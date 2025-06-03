# docs/conf.py
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
# If your project is in a subdirectory relative to `docs`,
# adjust the path accordingly.
# For example, if your project root is one level up from docs/:
sys.path.insert(0, os.path.abspath('..'))
# If your project root is two levels up (e.g., docs/source/):
# sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'My Awesome Project'
copyright = '2024, Your Name'
author = 'Your Name'
release = '0.1.0' # The short X.Y version
version = '0.1.0' # The full version, including alpha/beta/rc tags

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',     # For documenting Python code automatically
    'sphinx.ext.napoleon',    # For Google style or NumPy style docstrings
    'sphinx.ext.todo',        # For including todo notes in the documentation
    'sphinx.ext.viewcode',    # For linking to source code from docs
    'sphinx.ext.intersphinx', # For linking to other Sphinx docs (e.g., Python, NumPy)
    'sphinx_rtd_theme',       # Theme
    'myst_parser',            # If you want to write documentation in Markdown
]

# Set the source file suffixes. For Markdown and reStructuredText.
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# The master toctree document.
root_doc = 'index' # If your main document is index.rst/index.md

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme' # Use the Read the Docs theme

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for sphinx.ext.todo ---------------------------------------------
todo_include_todos = True # Show todo notes in the documentation
