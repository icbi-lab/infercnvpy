---
jupyter:
  jupytext:
    formats: md,ipynb
    notebook_metadata_filter: -kernelspec
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.0.rc1
---

# Infer CNV on lung cancer dataset

<!-- #raw raw_mimetype="text/restructuredtext" -->
In this tutorial, we demonstrate how infercnvpy can be used to derive Copy Number Variation (CNV)
from a single-cell RNA-seq dataset and to distinguish between tumor and normal cells. 

For this tutorial, we use the dataset by :cite:`Maynard2020` generated on the Smart-seq2 platform. 
The original dataset contains about 20,000 cells. Here, we use a downsampled version with 3,000 cells which is available through :func:`infercnvpy.datasets.maynard2020_3k`.

.. warning::
    We consider this method still experimental, but decided to already put it on GitHub because it might be useful. 
    Treat the results with care. In particular, no validation with ground-truth data has been performed. 
<!-- #endraw -->

```python
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt

sc.settings.set_figure_params(figsize=(5, 5))
```

```python
sc.logging.print_header()
```

## Loading the example dataset

<!-- #raw raw_mimetype="text/restructuredtext" -->
.. note::

    **Preprocessing data**
    
    Low-quality cells should already be filtered out and the input data 
    must be normalized and log-transformed. For more information, see
    :ref:`input-data`. 
    
    Also, the genomic positions need to be stored in `adata.var`. The 
    columns `chromosome`, `start`, and `end` hold the chromosome and 
    the start and end positions on that chromosome for each gene, 
    respectively. 
    
    Infercnvpy provides the :func:`infercnvpy.io.genomic_position_from_gtf` function
    to read these information from a GTF file and add them to `adata.var`. 
    
The example dataset is already appropriately preprocessed. 
<!-- #endraw -->

```python
adata = cnv.datasets.maynard2020_3k()
adata.var.loc[:, ["ensg", "chromosome", "start", "end"]].head()
```

Let's first inspect the UMAP plot based on the transcriptomics data:

```python
sc.pl.umap(adata, color="cell_type")
```

## Running infercnv

<!-- #raw raw_mimetype="text/restructuredtext" -->
Let's now run :func:`infercnvpy.tl.infercnv`. Essentially, this method sorts genes
by chromosome and genomic position and compares the average gene expression over genomic
region to a reference. The original inferCNV method uses a window size of 100, 
but larger window sizes can make sense, depending on the number of 
genes in your dataset. 

:func:`~infercnvpy.tl.infercnv` adds a `cell x genomic_region` matrix to 
`adata.obsm["X_cnv"]`. 

For more information about the method check out :ref:`infercnv-method`. 

.. note::

    **Choosing reference cells**
    
    The most common use-case is to compare tumor against normal cells. If you have 
    prior information about which cells are normal (e.g. from cell-type annotations
    based on transcriptomics data), it is recommended to provide this information 
    to :func:`~infercnvpy.tl.infercnv`. 
    
    The more different cell-types you 
    can provide, the better. Some cell-types physiologicallly over-express
    certain genomic regions (e.g. plasma cells highly express Immunoglobulin genes
    which are genomically adjacent). If you provide multiple cell-types, 
    only regions are considered being subject to CNV that are different from *all*
    provided cell-types. 
    
    If you don't provide any reference, the mean of all cells is used instead,
    which may work well on datasets that contain enough tumor and normal cells. 
<!-- #endraw -->

```python
# We provide all immune cell types as "normal cells".
cnv.tl.infercnv(
    adata,
    reference_key="cell_type",
    reference_cat=[
        "B cell",
        "Macrophage",
        "Mast cell",
        "Monocyte",
        "NK cell",
        "Plasma cell",
        "T cell CD4",
        "T cell CD8",
        "T cell regulatory",
        "mDC",
        "pDC",
    ],
    window_size=250,
)
```

Now, we can plot smoothed gene expression by cell-type and chromosome. 
We can observe that the Epithelial cell cluster, which consists largely of tumor cells, appears
to be subject to copy number variation. 

```python
cnv.pl.chromosome_heatmap(adata, groupby="cell_type")
```

## Clustering by CNV profiles and identifying tumor cells

<!-- #raw raw_mimetype="text/restructuredtext" -->
To cluster and annotate cells, `infercnvpy` mirrors the `scanpy` workflow. 
The following functions work exactely as their `scanpy` counterpart, except that 
they use the CNV profile matrix as input. Using these functions, we can perform
graph-based clustering and generate a UMAP plot based on the CNV profiles. 
Based on these clusters, we can annotate tumor and normal cells. 

.. module:: infercnvpy
   :noindex:

.. autosummary::
       
   infercnvpy.tl.pca
   infercnvpy.pp.neighbors
   infercnvpy.tl.leiden
   infercnvpy.tl.umap
   infercnvpy.pl.umap
<!-- #endraw -->

```python
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
```

After running leiden clustering, we can plot the chromosome heatmap 
by CNV clusters. We can observe that, as opposted to the clusters 
at the bottom, the clusters at the top have essentially no differentially expressed genomic regions. 
The differentially expressed regions are likely due to copy number variation and the respective 
clusters likely represent tumor cells. 

```python
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True)
```

### UMAP plot of CNV profiles

<!-- #raw raw_mimetype="text/restructuredtext" -->
We can visualize the same clusters as a UMAP plot. Additionally, 
:func:`infercnvpy.tl.cnv_score` computes a summary score that quantifies the amount of copy
number variation per cluster. It is simply defined as the
mean of the absolute values of the CNV matrix for each cluster. 
<!-- #endraw -->

```python
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)
```

The UMAP plot consists of a large blob of "normal" cells and several smaller clusters
with distinct CNV profiles. Except for cluster "12", which consists of ciliated cells, 
the isolated clusters are all epithelial cells. These are likely tumor cells and each 
cluster represents an individual sub-clone.

```python
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    adata,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata, color="cell_type", ax=ax3)
```

We can also visualize the CNV score and clusters on the transcriptomics-based UMAP plot. 
Again, we can see that there are subclusters of epithelial cells that belong
to a distinct CNV cluster, and that these clusters tend to have the 
highest CNV score. 

```python
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    2, 2, figsize=(12, 11), gridspec_kw=dict(wspace=0.5)
)
ax4.axis("off")
sc.pl.umap(adata, color="cnv_leiden", ax=ax1, show=False)
sc.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
sc.pl.umap(adata, color="cell_type", ax=ax3)
```

### Classifying tumor cells

Based on these observations, we can now assign cell to either "tumor" or "normal". 
To this end, we add a new column `cnv_status` to `adata.obs`. 

```python
adata.obs["cnv_status"] = "normal"
adata.obs.loc[
    adata.obs["cnv_leiden"].isin(["10", "13", "15", "5", "11", "16", "12"]), "cnv_status"
] = "tumor"
```

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw=dict(wspace=0.5))
cnv.pl.umap(adata, color="cnv_status", ax=ax1, show=False)
sc.pl.umap(adata, color="cnv_status", ax=ax2)
```

Now, we can plot the CNV heatmap for tumor and normal cells separately: 

```python
cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "tumor", :])
```

```python
cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "normal", :])
```
