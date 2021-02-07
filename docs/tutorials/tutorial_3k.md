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

For this tutorial, we use the dataset by :cite:t:`Maynard2020` generated on the Smart-seq2 platform. 
The original dataset contains about 20,000 cells. Here, we use a downsampled version with 3,000 cells which is available through :func:`infercnvpy.datasets.maynard2020_3k`.

.. warning::
    This tutorial is still work-in-progress!
<!-- #endraw -->

```python
import scanpy as sc
import infercnvpy as cnv
```

```python
sc.logging.print_header()
```

```python
# adata = cnv.datasets.maynard2020_3k()
```

```python

```
