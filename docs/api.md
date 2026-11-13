# API

Import infercnvpy together with scanpy as

```python
import scanpy as sc
import infercnvpy as cnv
```

For consistency, the infercnvpy API tries to follow the {doc}`scanpy API <scanpy:api/index>` as closely as possible.

(api-io)=

## Input/Output: `io`

```{eval-rst}
.. module:: infercnvpy.io
.. currentmodule:: infercnvpy

.. autosummary::
    :toctree: generated

    io.genomic_position_from_biomart
    io.genomic_position_from_gtf
    io.read_scevan
```

## Preprocessing: `pp`

```{eval-rst}
.. module:: infercnvpy.pp
.. currentmodule:: infercnvpy

.. autosummary::
    :toctree: generated

    pp.neighbors
```

## Tools: `tl`

Tools add an interpretable annotation to the {class}`~anndata.AnnData` object,
which usually can be visualized by a corresponding plotting function.

The tools for embeddings and clustering mirror the scanpy API.
However, while the scanpy tools operate on transcriptomics data,
the infercnvpy equivalents operate on CNV data.

```{eval-rst}
.. module:: infercnvpy.tl
.. currentmodule:: infercnvpy
```

### InferCNV

```{eval-rst}
.. autosummary::
    :toctree: generated

    tl.infercnv
    tl.copykat
```

### CNV scores

```{eval-rst}
.. autosummary::
    :toctree: generated

    tl.cnv_score
    tl.ithcna
    tl.ithgex
```

### Embeddings

```{eval-rst}
.. autosummary::
    :toctree: generated

    tl.pca
    tl.umap
    tl.tsne
```

### Clustering

```{eval-rst}
.. autosummary::
    :toctree: generated

    tl.leiden
```

## Plotting: `pl`

```{eval-rst}
.. module:: infercnvpy.pl
.. currentmodule:: infercnvpy
```

### InferCNV

```{eval-rst}
.. autosummary::
    :toctree: generated

    pl.chromosome_heatmap
    pl.chromosome_heatmap_summary
```

### Embeddings

```{eval-rst}
.. autosummary::
    :toctree: generated

    pl.umap
    pl.tsne
```

## Datasets: `datasets`

```{eval-rst}
.. module:: infercnvpy.datasets
.. currentmodule:: infercnvpy

.. autosummary::
    :toctree: generated

    datasets.maynard2020_3k
    datasets.oligodendroglioma
```
