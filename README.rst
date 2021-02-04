infercnvpy: Scanpy plugin to infer copy number variation (CNV) from single-cell transcriptomics data
====================================================================================================
|tests| |docs| |pypi| |black|

.. |tests| image:: https://github.com/grst/infercnvpy/workflows/tests/badge.svg
    :target: https://github.com/grst/infercnvpy/actions?query=workflow%3Atests
    :alt: Build Status

.. |docs| image::  https://github.com/grst/infercnvpy/workflows/docs/badge.svg
    :target: https://grst.github.io/infercnvpy
    :alt: Documentation Status

.. |pypi| image:: https://img.shields.io/pypi/v/infercnvpy?logo=PyPI
    :target: https://pypi.org/project/infercnvpy/
    :alt: PyPI

.. .. |bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
..      :target: http://bioconda.github.io/recipes/infercnvpy/README.html
..      :alt: Bioconda

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: The uncompromising python formatter

Infercnv is a scalable python library to infer copy number variation (CNV) events
from single cell transcriptomics data. It is heavliy inspired by `InferCNV <https://github.com/broadinstitute/inferCNV/wiki>`_,
but plays nicely with `scanpy <https://scanpy.readthedocs.io/en/stable/index.html>`_ and is much more scalable.

.. image:: img/infercnv_heatmap.png
    :align: center
    :alt: The main result of infercnv

Getting started
^^^^^^^^^^^^^^^

**WARNING: This package is still experimental. The results have not been validated,
except in that they look similar, but not identical, to the results of InferCNV.**

Please refer to the `documentation <https://grst.github.io/infercnvpy>`_. In particular, the

- `Tutorial <https://grst.github.io/infercnvpy/tutorials/tutorial_3k.html>`_, and the
- `API documentation <https://grst.github.io/infercnvpy/api.html>`_.


Installation
^^^^^^^^^^^^
You need to have Python 3.8 or newer installed on your system. If you don't have
Python installed, we recommend installing `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.

There are several alternative options to install infercnvpy:

1) Install the latest release of `infercnvpy` from `PyPI <https://pypi.org/project/infercnvpy/>`_:

.. code-block::

    pip install infercnvpy


.. 2) Get it from `Bioconda <http://bioconda.github.io/recipes/infercnvpy/README.html>`_:

.. .. code-block::

..     conda install -c conda-forge -c bioconda infercnvpy


2) Install the latest development version:

.. code-block::

    pip install git+https://github.com/grst/infercnvpy.git@master


.. 4) Run it in a container using `Docker <https://www.docker.com/>`_ or `Podman <https://podman.io/>`_:

.. .. code-block::

..     docker pull quay.io/biocontainers/infercnvpy:<tag>

.. where `tag` is one of `these tags <https://quay.io/repository/biocontainers/infercnvpy?tab=tags>`_.


Release notes
^^^^^^^^^^^^^
See the `release section <https://github.com/grst/infercnvpy/releases>`_.

Contact
^^^^^^^
Please use the `issue tracker <https://github.com/grst/infercnvpy/issues>`_.

Citation
^^^^^^^^
n/a
