"""Python library to infer copy number variation (CNV) from single-cell RNA-seq data"""

from ._metadata import __version__, __author__, __email__, within_flit

if not within_flit():
    from . import io, pl, pp, tl, datasets
