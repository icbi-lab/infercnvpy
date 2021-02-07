---
jupyter:
  jupytext:
    notebook_metadata_filter: -kernelspec
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.0.rc1
---

# Reproduce the heatmap from inferCNV

This document demonstrates to reproduce how the [example heatmap](https://github.com/broadinstitute/inferCNV/wiki#demo-example-figure) from the original
R inverCNV implementation. It is based on a small, 183-cell example dataset of malignant and non-malignant cells from Oligodendroglioma derived from [Tirosh et al. (2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5465819/). 

```python
import infercnvpy as cnv
import scanpy as sc
```

## Prepare and inspect dataset
The example dataset is available in the `datasets` module. It is already TPM-normalized, but not log-transformed. 

```python
adata = cnv.datasets.oligodendroglioma()
sc.pp.log1p(adata)
```

It also already has the genomic positions annotated in`adata.var`: 

```python
adata.var.head()
```

It contains four types of malignant cells, and two clusters of non-malignant cells. 

```python
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color="cell_type")
```

## Run infercnvpy

In this case we know which cells are non-malignant. For best results, it is recommended to use
the non-malignant cells as a background. We can provide this information using `reference_key` and `reference_cat`. 

In order to reproduce the results as exactely as possible, we use a `window_size` of 100 and a `step` of 1. 

```python
%%time
cnv.tl.infercnv(
    adata, 
    reference_key="cell_type", 
    reference_cat=["Oligodendrocytes (non-malignant)", "Microglia/Macrophage"], 
    window_size=100,
    step=1,
)
```

```python
%%time
cnv.pl.chromosome_heatmap(adata, groupby="cell_type", dendrogram=True)
```

Note that running the same analysis in R (`invercnv v1.6.0` from Bioconductor) takes about 1:30 min. 
