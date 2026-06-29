# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

project = 'PyMERK'
copyright = '2026, Pierre Beaujean'
author = 'Pierre Beaujean'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.append(os.path.abspath('extensions'))

extensions = [
    'sphinx.ext.doctest',
    'sphinx.ext.napoleon',
    'sphinx.ext.githubpages',
    'sphinx.ext.mathjax',
    'sphinxcontrib.mermaid',
    # custom extension for this project
    'pym_keywords',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Napoleon settings
napoleon_google_docstring = True
napoleon_include_init_with_doc = True

# Mermaid
mermaid_d3_zoom = True
mermaid_height = "1000px"

# -- Options for HTML output -------------------------------------------------
html_theme = "shibuya"

html_theme_options = {
  "github_url": "https://github.com/pierre-24/pymerk"
}

html_context = {
    "source_type": "github",
    "source_user": "pierre-24",
    "source_repo": "pymerk",
    "source_edit_template": "https://github.com/pierre-24/pymerk/blob/main/docs/{0}",
}

html_sidebars = {
  "**": [
    "sidebars/localtoc.html",
    "sidebars/repo-stats.html",
    "sidebars/edit-this-page.html",
  ]
}