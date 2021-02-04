"""Example datasets"""
import importlib.resources as pkg_resources
from . import data
import scanpy as sc
from anndata import AnnData


def oligodendroglioma() -> AnnData:
    """The original inferCNV example dataset.

    Derived from :cite:t:`Tirosh2016`.
    """
    with pkg_resources.path(data, "oligodendroglioma.h5ad") as p:
        return sc.read_h5ad(p)


def maynard2020_3k() -> AnnData:
    raise NotImplementedError
