"""Python library to infer copy number variation (CNV) from single-cell RNA-seq data"""

from ._metadata import __author__, __email__, __version__, within_flit

if not within_flit():
    from . import datasets, io, pl, pp, tl
