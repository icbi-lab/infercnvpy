# Configuration file for the Sphinx documentation builder.
# `scverse_doc` sets the theme and shared defaults, see https://scverse-doc.readthedocs.io
# For the full list of options, see https://www.sphinx-doc.org/page/usage/configuration.html

from datetime import datetime
from importlib.metadata import metadata

from scverse_doc.registry import intersphinx

# NOTE: If you installed your project in editable mode, this might be stale.
#       If this is the case, reinstall it to refresh the metadata
info = metadata("infercnvpy")
project = info["Name"]
author = info["Author"]
copyright = f"{datetime.now():%Y}, {author}."
version = release = info["Version"]
urls = dict(pu.split(", ") for pu in info.get_all("Project-URL"))
source_repository = urls["Source"]

extensions = [
    "scverse_doc",
    "scverse_misc.sphinx_ext",
    "sphinxcontrib.bibtex",
    "sphinx_autodoc_typehints",
    "IPython.sphinxext.ipython_console_highlighting",
]

bibtex_bibfiles = ["references.bib"]
typehints_defaults = "braces"
always_use_bars_union = True  # use `|` instead of `Union` in types even when building with Python ≤3.14

intersphinx_mapping = intersphinx() | {"pooch": ("https://www.fatiando.org/pooch/latest/", None)}

html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

xref_overrides = {
    # If building the documentation fails because a type annotation links to a private path
    # (e.g. `igraph._igraph.GraphBase`), map it to the documented name here.
    # For a missing link that can’t be fixed this way, add it to `nitpick_ignore` instead.
    #     "igraph._igraph.GraphBase": "igraph.GraphBase",
}

# lists the adata fields a function writes, rendered as a definition list
napoleon_custom_sections = ["Sets the following fields"]
