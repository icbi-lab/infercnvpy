"""Python library to infer copy number variation (CNV) from single-cell RNA-seq data."""

from importlib.metadata import version

from . import datasets, io, pl, pp, tl

__all__ = ["datasets", "io", "pl", "pp", "tl"]
__version__ = version("infercnvpy")
