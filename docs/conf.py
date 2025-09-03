import sphinx_autosummary_accessors
import shnitsel.xarray

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'shnitsel-tools'
copyright = '2025, Theodor Everley Röhrkasten, Carolin Müller'
author = 'Theodor Everley Röhrkasten, Carolin Müller'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',               # Reads docstrings from sourcecode
    'sphinx.ext.napoleon',              # Converts numpy-formatted docstrings to reST
    'sphinx_autosummary_accessors',     
    'sphinx.ext.intersphinx',           # Produces links to external Sphinx-based docs
    'autoapi.extension'                 # Finds all things to document automatically
]

templates_path = ['_templates', sphinx_autosummary_accessors.templates_path]
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'test_autosummary_accessors.rst']

add_module_names = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# -- Options for autodoc extension -------------------------------------------

autodoc_typehints = 'description'
# napoleon_use_param = False
autodoc_type_aliases = {
    'ArrayLike': 'numpy.typing.ArrayLike',
}

# -- Options for intersphinx extension ---------------------------------------
intersphinx_mapping = {
    'xarray': ('https://docs.xarray.dev/en/stable/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'python': ('https://docs.python.org/3', None),
}

# -- Options for autoapi extension -------------------------------------------
autoapi_dirs = ['../shnitsel']
autoapi_root = 'api/'
# autoapi_keep_files = True