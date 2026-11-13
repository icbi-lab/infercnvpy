(large-data)=

# Large and lazy datasets

{func}`infercnvpy.tl.infercnv` processes the cells in blocks and never needs the whole expression matrix in memory
as a dense array. The expression matrix can be

- an in-memory `numpy` array or `scipy` sparse matrix,
- a [dask](https://www.dask.org/) array, also with sparse `scipy` blocks, e.g. a matrix loaded lazily with
  {func}`anndata.experimental.read_elem_lazy` or {func}`anndata.experimental.read_lazy`, or a table from
  [spatialdata-io](https://spatialdata.scverse.org/projects/io/),
- the matrix of a backed AnnData (`anndata.read_h5ad(..., backed="r")`).

Install the optional dependency with `pip install "infercnvpy[dask]"` to use dask arrays.

## How the data is read

infercnv reads the data up to three times, one block of cells at a time:

1. For a dask array or a backed matrix, the reference expression (the mean of all cells, or of every reference
   category) is computed in one pass over the data. In-memory data needs no extra pass.
2. The noise band of `dynamic_threshold` is computed from the blocks that hold reference cells, or from all blocks
   if there are no reference cells or `noise_reference="all"`. This pass is skipped when all cells fit in one block.
3. The CNV matrix is computed block by block.

In memory and backed, a block has `chunksize` cells. A dask array keeps its own blocks, and `chunksize` is ignored.
The block size does not change the result, except that the reference means of lazy data are sums over blocks, so
their last digits can differ from the in-memory means. Peak memory is the result plus about `chunksize` x windows
float32 values for the block being processed (and, for lazy data, the block that was read).

## Threads and dask schedulers

By default, the blocks are read one after another on the main thread, with the synchronous scheduler for a dask
array, and every block is processed with `n_jobs` threads. The result, `adata.obsm["X_cnv"]`, is an in-memory
sparse matrix with 8 bytes per value that is not removed as noise. With 5,000 windows of which 10% are kept, that is
about 4 GB per million cells.

With `lazy=True`, the result for a dask array stays a dask array of CSR blocks, computed with
{meth}`dask.array.Array.map_blocks`. The noise band is computed first, with one dask computation. When a
{class}`distributed.Client` is active, `lazy=True` is the default. In this mode the dask scheduler decides how
many blocks run at the same time, each in a single thread, and `n_jobs` is ignored. infercnvpy never starts numba
threads inside a dask worker thread, because numba's default threading layer on macOS aborts the process when it is
used from several threads at once.

:::{note}
infercnv expects log-normalised expression (see {ref}`input-data`). For a dask array, {func}`scanpy.pp.normalize_total`
and {func}`scanpy.pp.log1p` stay lazy. A backed AnnData cannot be normalised in place by scanpy, so normalise the data
before you write the file.
:::

## Example: a lazily loaded table

Spatial technologies often produce tables with a lazy dask `X` of CSR blocks that are indexed by Ensembl ID. Next to
the genes, such a table can contain control probes or codewords without a genomic position, which you remove before
normalisation.

```python
import anndata as ad
import scanpy as sc
import zarr
from anndata.experimental import read_elem_lazy

import infercnvpy as cnv

# `obs` and `var` in memory, `X` as a dask array (or use the table of a SpatialData object)
f = zarr.open("table.zarr", mode="r")
table = ad.AnnData(X=read_elem_lazy(f["X"]), obs=ad.io.read_elem(f["obs"]), var=ad.io.read_elem(f["var"]))

# genomic positions: the table is indexed by Ensembl ID, so match the GTF's `gene_id`
cnv.io.genomic_position_from_gtf("gencode.v44.annotation.gtf", table, gtf_gene_id="gene_id")

# keep gene features only, and normalise lazily
genes = table[:, (table.var["feature_type"] == "Gene Expression").to_numpy()].copy()
sc.pp.normalize_total(genes, target_sum=1e4)
sc.pp.log1p(genes)

# `cell_type` and the categories stand for an annotation of normal cells in `obs`
cnv.tl.infercnv(genes, reference_key="cell_type", reference_cat=["Fibroblast", "Immune"])
```

Genes on chromosomes without a position (and on `chrX` and `chrY` by default) are skipped, and a chromosome with
fewer genes than `window_size` becomes a single window. Targeted panels with few genes per chromosome carry much less
CNV signal than whole-transcriptome data, so consider a smaller `window_size`.

## Plotting

{func}`infercnvpy.pl.chromosome_heatmap` averages the cells that fall into each pixel row of the figure, so it reads
the CNV matrix once and draws an image of about a thousand rows, whatever the number of cells. It needs
`adata.obsm["X_cnv"]` in memory: compute a lazy result first, e.g.
`adata.obsm["X_cnv"] = adata.obsm["X_cnv"].compute()`.
