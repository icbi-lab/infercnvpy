# Installation

infercnvpy needs Python 3.13 or newer.
We recommend [uv][], which installs Python for you and keeps each project's dependencies in their own virtual environment.

## PyPI

Install the latest release from [PyPI][]:

:::::{tab-set}
::::{tab-item} uv
:sync: uv

```bash
uv add infercnvpy
```

::::

::::{tab-item} Pip
:sync: pip

```bash
pip install infercnvpy
```

::::
:::::

## Development version

Install the latest development version from GitHub:

:::::{tab-set}
::::{tab-item} uv
:sync: uv

```bash
uv add git+https://github.com/scverse/infercnvpy.git
```

::::

::::{tab-item} Pip
:sync: pip

```bash
pip install git+https://github.com/scverse/infercnvpy.git
```

::::
:::::

## Optional dependencies

- `infercnvpy[gtf]` installs [gtfparse][] for {func}`infercnvpy.io.genomic_position_from_gtf`.
- `infercnvpy[copykat]` installs [rpy2][] for {func}`infercnvpy.tl.copykat`.
  It also needs a working R installation with the [copykat][] package.
  Usually, if `R` is in your `PATH`, rpy2 finds your R installation.
  If importing rpy2 fails, set the `R_HOME` environment variable before importing infercnvpy:

  ```python
  import os

  os.environ["R_HOME"] = "/usr/lib/R"
  import infercnvpy
  ```

To work on infercnvpy itself, see {doc}`contributing`.

[uv]: https://docs.astral.sh/uv/
[pypi]: https://pypi.org/project/infercnvpy
[gtfparse]: https://github.com/openvax/gtfparse
[rpy2]: https://rpy2.github.io/
[copykat]: https://github.com/navinlabcode/copykat#step-1-installation
