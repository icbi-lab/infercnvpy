# Quickstart

infercnvpy works on {class}`~anndata.AnnData` objects with log-normalized expression and the genomic position of each gene in `adata.var`.
{mod}`infercnvpy.tl` infers copy number variation and clusters cells by it, {mod}`infercnvpy.pp` builds the neighborhood graph and {mod}`infercnvpy.pl` plots the result.

```python
import scanpy as sc

import infercnvpy as cnv

adata = cnv.datasets.oligodendroglioma()

cnv.tl.infercnv(
    adata,
    reference_key="cell_type",
    reference_cat=["Oligodendrocytes (non-malignant)", "Microglia/Macrophage"],
)
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.cnv_score(adata)

cnv.pl.chromosome_heatmap(adata, groupby="cell_type")
```

If `adata.var` has no `chromosome`, `start` and `end` columns yet, add them with {func}`infercnvpy.io.genomic_position_from_gtf` or {func}`infercnvpy.io.genomic_position_from_biomart`.
{doc}`infercnv` describes the method and the input it needs, the {doc}`tutorials <tutorials>` go into more detail and the {doc}`API reference <api>` lists every function.
