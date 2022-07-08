# infercnvpy

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]
[![PyPI][badge-pypi]][link-pypi]

[badge-tests]: https://img.shields.io/github/workflow/status/icbi-lab/infercnvpy/Test/main
[link-tests]: https://github.com/icbi-lab/infercnvpy/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/infercnvpy
[badge-pypi]: https://img.shields.io/pypi/v/infercnvpy?logo=PyPI
[link-pypi]: https://pypi.org/project/infercnvpy/

Infercnv is a scalable python library to infer copy number variation (CNV) events from single cell transcriptomics data. It is heavliy inspired by InferCNV, but plays nicely with scanpy and is much more scalable.

![The main result of infercnv](img/infercnv_heatmap.png)

**WARNING**:

**This package is still experimental. The results have not been validated,
except in that they look similar, but not identical, to the results of InferCNV.**

**We are happy about feedback and welcome contributions!**

## Getting started

Please refer to the [documentation][link-docs]. In particular, the

-   [API documentation][link-api].

## Installation

You need to have Python 3.8 or newer installed on your system. If you don't have
Python installed, we recommend installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

There are several alternative options to install infercnvpy:

1. Install the latest release of `infercnvpy` from `PyPI <https://pypi.org/project/infercnvpy/>`\_:

```bash
pip install infercnvpy
```

2. Install the latest development version:

```bash
pip install git+https://github.com/icbi-lab/infercnvpy.git@main
```

## Release notes

See the [changelog][changelog].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][scverse-discourse].
If you found a bug, please use the [issue tracker][issue-tracker].

## Citation

> t.b.a

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/icbi-lab/infercnvpy/issues
[changelog]: https://infercnvpy.readthedocs.io/latest/changelog.html
[link-docs]: https://infercnvpy.readthedocs.io
[link-api]: https://infercnvpy.readthedocs.io/latest/api.html
